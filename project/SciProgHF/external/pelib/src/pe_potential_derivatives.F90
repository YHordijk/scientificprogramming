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
module pe_potential_derivatives

    use pe_precision

    implicit none

    private

    public :: electron_fields
    public :: nuclear_fields
    public :: multipole_fields

contains

!------------------------------------------------------------------------------

subroutine electron_fields(Fels, denmats)

    use pe_mpi
    use pe_variables
    use pe_integral_interfaces
    use pe_blas_interfaces, only: nrm2, dot

    real(dp), dimension(:,:), intent(out) :: Fels
    real(dp), dimension(:), intent(in) :: denmats

    logical :: skip
    integer :: site
    integer :: i, j, k, l, m
    real(dp) :: factor
    real(dp), dimension(3) :: Rms
    real(dp), dimension(:,:,:), allocatable :: Fel_ints

    allocate(Fel_ints(nnbas,3,1))

    i = 0
    do site = site_start, site_finish
        if (zeroalphas(site)) cycle
        call Tk_integrals('potential_derivative', Rs(:,site), Fel_ints)
        do j = 1, 3
            do k = 1, ndens
                l = (k - 1) * nnbas + 1
                m = k * nnbas
                Fels(i+j,k) = dot(denmats(l:m), Fel_ints(:,j,1))
            end do
        end do
        if (pe_core_damp) then
            call get_damping_factor(Rs(:,site), P11s(:,site), factor)
            do k = 1, ndens
                Fels(i+1:i+3,k) = factor * Fels(i+1:i+3,k)
            end do
        endif
        i = i + 3
    end do

#if defined(VAR_MPI)
    if (myid == master .and. nprocs > 1) then
        displs(0) = 0
        do i = 1, nprocs
            displs(i) = displs(i-1) + poldists(i-1)
        end do
        do i = 1, ndens
            call mpi_gatherv(mpi_in_place, dummy_int, rmpi, Fels(1,i), poldists, displs, rmpi, master, comm, ierr)
        end do
    else if (myid /= master) then
        do i = 1, ndens
            call mpi_gatherv(Fels(1,i), poldists(myid), rmpi, dummy_real, poldists, displs, rmpi, master, comm, ierr)
        end do
    end if
#endif

    deallocate(Fel_ints)

end subroutine electron_fields

!------------------------------------------------------------------------------

subroutine nuclear_fields(Fnucs)

    use pe_mpi
    use pe_variables
    use pe_utils, only: openfile
    use pe_multipole_tools, only: Tk_tensor
    use pe_blas_interfaces, only: nrm2

    real(dp), dimension(:), intent(out) :: Fnucs

    logical :: lexist, skip
    integer :: lu, site
    integer :: i, j, k
    real(dp) :: factor
    real(dp), dimension(3) :: Rms, Tms

    if (myid == master) then
        inquire(file='pe_nuclear_field.bin', exist=lexist)
    end if

#if defined(VAR_MPI)
    if (nprocs > 1) then
        call mpi_bcast(lexist, 1, lmpi, master, comm, ierr)
    end if
#endif

    if (lexist .and. ((scf_cycle > 1) .or. pe_restart)) then
        if (myid == master) then
            call openfile('pe_nuclear_field.bin', lu, 'old', 'unformatted')
            rewind(lu)
            read(lu) Fnucs
            close(lu)
        end if
    else
        Fnucs = 0.0_dp
        i = 0
        do site = site_start, site_finish
            if (zeroalphas(site)) cycle
            do j = 1, nnucs
                Rms = Rs(:,site) - Rm(:,j)
                call Tk_tensor(Tms, Rms)
                do k = 1, 3
                    Fnucs(i+k) = Fnucs(i+k) - Zm(1,j) * Tms(k)
                end do
            end do
            if (pe_core_damp) then
                call get_damping_factor(Rs(:,site), P11s(:,site), factor)
                Fnucs(i+1:i+3) = factor * Fnucs(i+1:i+3)
            end if
            i = i + 3
        end do
#if defined(VAR_MPI)
        if (myid == master .and. nprocs > 1) then
            displs(0) = 0
            do i = 1, nprocs
                displs(i) = displs(i-1) + poldists(i-1)
            end do
            call mpi_gatherv(mpi_in_place, dummy_int, rmpi, Fnucs(1), poldists, displs, rmpi, master, comm, ierr)
        else if (myid /= master) then
            call mpi_gatherv(Fnucs(1), poldists(myid), rmpi, dummy_real, poldists, displs, rmpi, master, comm, ierr)
        end if
#endif
        if (myid == master) then
            call openfile('pe_nuclear_field.bin', lu, 'unknown', 'unformatted')
            rewind(lu)
            write(lu) Fnucs
            close(lu)
        end if
    end if

end subroutine nuclear_fields

!------------------------------------------------------------------------------

subroutine multipole_fields(Fmuls)

    use pe_mpi
    use pe_variables
    use pe_utils, only: openfile
    use pe_multipole_tools, only: multipole_derivative, multipole_derivative_damped

    real(dp), dimension(:), intent(out) :: Fmuls

    logical :: exclude, lexist
    integer :: lu
    integer :: i, j, k, l
    real(dp) :: alpha_i, alpha_j
    real(dp), dimension(3) :: Rji, Fi
    real(dp), parameter :: i3 = 1.0_dp / 3.0_dp

    if (myid == master) then
        inquire(file='pe_multipole_field.bin', exist=lexist)
    end if

#if defined(VAR_MPI)
    if (nprocs > 1) then
        call mpi_bcast(lexist, 1, lmpi, master, comm, ierr)
    end if
#endif
    if (lexist .and. ((scf_cycle > 1) .or. pe_restart)) then
        if (myid == master) then
            call openfile('pe_multipole_field.bin', lu, 'old', 'unformatted')
            rewind(lu)
            read(lu) Fmuls
            close(lu)
        end if
    else
        Fmuls = 0.0_dp
        l = 1
        do i = 1, nsites
            if (zeroalphas(i)) cycle
            if (pe_mul_damp) then
                alpha_i = (P11s(1,i) + P11s(4,i) + P11s(6,i)) * i3
            end if
            do j = site_start, site_finish
                if (i == j) then
                    cycle
                end if
                exclude = .false.
                do k = 1, lexlst
                    if (exclists(k,i) == exclists(1,j)) then
                        exclude = .true.
                        exit
                    end if
                end do
                if (exclude) cycle
                if (pe_mul_damp) then
                    alpha_j = (P11s(1,j) + P11s(4,j) + P11s(6,j)) * i3
                end if
! TODO: cutoff???
                Fi = 0.0_dp
                Rji = Rs(:,i) - Rs(:,j)
                if (lmul(0)) then
                    if (maxval(abs(M0s(:,j))) >= zero) then
                        if (pe_mul_damp) then
                            call multipole_derivative_damped(Fi, Rji, M0s(:,j), alpha_i, alpha_j, mul_damp)
                        else
                            call multipole_derivative(Fi, Rji, M0s(:,j))
                        end if
                    end if
                end if
                if (lmul(1)) then
                    if (maxval(abs(M1s(:,j))) >= zero) then
                        if (pe_mul_damp) then
                            call multipole_derivative_damped(Fi, Rji, M1s(:,j), alpha_i, alpha_j, mul_damp)
                        else
                            call multipole_derivative(Fi, Rji, M1s(:,j))
                        end if
                    end if
                end if
                if (lmul(2)) then
                    if (maxval(abs(M2s(:,j))) >= zero) then
                        if (pe_mul_damp) then
                            call multipole_derivative_damped(Fi, Rji, M2s(:,j), alpha_i, alpha_j, mul_damp)
                        else
                            call multipole_derivative(Fi, Rji, M2s(:,j))
                        end if
                    end if
                end if
                if (lmul(3)) then
                    if (maxval(abs(M3s(:,j))) >= zero) then
                        if (pe_mul_damp) then
                            call multipole_derivative_damped(Fi, Rji, M3s(:,j), alpha_i, alpha_j, mul_damp)
                        else
                            call multipole_derivative(Fi, Rji, M3s(:,j))
                        end if
                    end if
                end if
                if (lmul(4)) then
                    if (maxval(abs(M4s(:,j))) >= zero) then
                        if (pe_mul_damp) then
                            stop 'ERROR: damping not implemented for fourth-order multipoles'
                        else
                            call multipole_derivative(Fi, Rji, M4s(:,j))
                        end if
                    end if
                end if
                if (lmul(5)) then
                    if (maxval(abs(M5s(:,j))) >= zero) then
                        if (pe_mul_damp) then
                            stop 'ERROR: damping not implemented for fifth-order multipoles'
                        else
                            call multipole_derivative(Fi, Rji, M5s(:,j))
                        end if
                    end if
                end if
                Fmuls(l:l+2) = Fmuls(l:l+2) + Fi
            end do
            l = l + 3
        end do
#if defined(VAR_MPI)
        if (myid == master .and. nprocs > 1) then
            call mpi_reduce(mpi_in_place, Fmuls(1), 3 * npols, rmpi, mpi_sum, master, comm, ierr)
        else if (myid /= master) then
            call mpi_reduce(Fmuls(1), dummy_real, 3 * npols, rmpi, mpi_sum, master, comm, ierr)
        end if
#endif
        if (myid == master) then
            call openfile('pe_multipole_field.bin', lu, 'unknown', 'unformatted')
            rewind(lu)
            write(lu) Fmuls
            close(lu)
        end if
     end if

end subroutine multipole_fields

!------------------------------------------------------------------------------

subroutine get_damping_factor(Ri, P11i, factor)

    use pe_variables
    use pe_multipole_tools, only: damping_coefficient
    use pe_blas_interfaces, only: nrm2

    real(dp), dimension(3), intent(in) :: Ri
    real(dp), dimension(6), intent(in) :: P11i
    real(dp), intent(out) :: factor

    integer :: j, jnuc
    real(dp) :: Rold, Rnew, alpha_i, alpha_m, v
    real(dp), dimension(3) :: Rmi
    real(dp), parameter :: i3 = 1.0_dp / 3.0_dp

    ! Thole damping
    ! JPC A 102 (1998) 2399 & Mol. Sim. 32 (2006) 471
    ! v = a * u , where
    ! a = 2.1304 (default)
    ! u = R / (alpha_i * alpha_j)**(1/6)

    ! locate the atom closest to Ri
    jnuc = - 1
    Rold = 1.0e9_dp
    do j = 1, nnucs
        Rmi = Ri - Rm(:,j)
        Rnew = nrm2(Rmi)
        if (Rnew <= Rold) then
            Rold = Rnew
            jnuc = j
        endif
    enddo

    if (jnuc == - 1) stop 'ERROR: Core damping failed because no atom was found'

    Rmi = Ri - Rm(:,jnuc)
    alpha_i = (P11i(1) + P11i(4) + P11i(6)) * i3
    alpha_m = (core_alphas(1,jnuc) + core_alphas(4,jnuc) + core_alphas(6,jnuc)) * i3

    call damping_coefficient(Rmi, alpha_i, alpha_m, core_damp, v)

    factor = 1.0_dp - (1.0_dp + v + 0.5_dp * v**2) * exp(- v)

end subroutine get_damping_factor

!------------------------------------------------------------------------------

end module pe_potential_derivatives
