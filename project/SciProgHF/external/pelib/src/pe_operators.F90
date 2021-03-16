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
module pe_operators

    use pe_precision

    implicit none

    private

    public :: get_fock, get_magnetic_gradient, get_molecular_gradient
    public :: get_effdipole

contains

!------------------------------------------------------------------------------

subroutine get_fock(denmats, fckmats, energies)

    use pe_mpi
    use pe_variables
    use pe_utils, only: openfile

    real(dp), dimension(:), intent(in) :: denmats
    real(dp), dimension(:), intent(out), optional :: fckmats
    real(dp), dimension(:), intent(out), optional :: energies

    integer :: i, lu
    logical :: es, pol

    es = .false.
    pol = .false.

    if (any(lmul) .and. .not. response) es = .true.
    if (pe_polar) pol = .true.

    if (fock .or. energy) then
        allocate(Epe(ndens))
        allocate(Ees(3,ndens))
        allocate(Epol(3,ndens))
        allocate(Elj(2))
        Epe = 0.0_dp
        Ees = 0.0_dp
        Epol = 0.0_dp
        Elj = 0.0_dp
    end if

    if (fock) then
        if (es) call electrostatic_operator(denmats, fckmats)
        if (pol) call polarization_operator(denmats, fckmats)
    else if (energy) then
        if (es) call electrostatic_operator(denmats)
        if (pol) call polarization_operator(denmats)
    else if (response) then
        if (pol) call polarization_operator(denmats, fckmats)
    end if

    if (fock .or. energy) then
        if (lvdw) call lj_energy
        if (myid == master) then
            do i = 1, ndens
                Epe(i) = sum(Ees(:,i)) + sum(Epol(:,i)) + sum(Elj(:))
            end do
            if (fock) energies = Epe
            if (energy) call print_info()
        end if
    end if

#if defined(VAR_MPI)
    if (fock .or. response) then
        if (myid == master .and. nprocs > 1) then
            call mpi_reduce(mpi_in_place, fckmats(1), ndens * nnbas, rmpi, mpi_sum, master, comm, ierr)
        else if (myid /= master) then
            call mpi_reduce(fckmats(1), dummy_real, ndens * nnbas, rmpi, mpi_sum, master, comm, ierr)
        end if
    end if
#endif

    if (fock) then
        if (myid == master) then
            call openfile('pe_fock.bin', lu, 'unknown', 'unformatted')
            rewind(lu)
            write(lu) fckmats
            close(lu)
        end if
    end if

    if (fock .or. energy) then
        deallocate(Epe)
        deallocate(Ees)
        deallocate(Epol)
        deallocate(Elj)
    end if

end subroutine get_fock

!------------------------------------------------------------------------------

subroutine print_info()

    use pe_variables
    use pe_utils, only: openfile
    use pe_blas_interfaces, only: nrm2

    integer :: i, j
    integer :: lu, nwarn
    real(dp) :: indnrm
    real(dp), dimension(3) :: indtot
    real(dp), dimension(:,:), allocatable :: Mkinds
    logical :: lexist

    write(luout, '(/4x,a)') '.--------------------------------------------------.'
    write(luout, '(4x,a)')  '| Final results from polarizable embedding library |'
    write(luout, '(4x,a/)') '`--------------------------------------------------´'

    if (ndens > 1) then
        write(luout, '(a)') 'INFO: only printing information for first&
                            & input density matrix'
    end if

    if (lpol(1,1)) then
        allocate(Mkinds(3*npols,1))
        inquire(file='pe_induced_moments.bin', exist=lexist)

        if (lexist) then
            call openfile('pe_induced_moments.bin', lu, 'old', 'unformatted')
            rewind(lu)
            read(lu) Mkinds
            close(lu)
        end if

        ! check induced dipoles
        nwarn = 0
        j = 1
        do i = 1, nsites
            if (zeroalphas(i)) cycle
            indnrm = nrm2(Mkinds(j:j+2,1))
            if (indnrm > 1.0_dp .and. (nwarn <= 10 .or. pe_verbose)) then
                write(luout, '(a,i6,a,f8.4)') 'WARNING: induced dipole at site ', i, ' has a magnitude of: ', indnrm
                nwarn = nwarn + 1
            else if (nwarn > 10) then
                write(luout, '(a)') 'WARNING: there are more than ten induced dipoles with magnitude > 1.0 a.u.'
                write(luout, '(a)') '         (use .VERBOSE to print all of them)'
                exit
            end if
            j = j + 3
        end do
        if (nwarn > 0) then
            write(luout, '(a)') 'WARNING: there are induced dipoles with magnitude larger than 1.0 a.u.'
            write(luout, '(a)') '         This can be an indication that there are problems with the structure or'
            write(luout, '(a)') '         something else that can cause over-polarization. Some over-polarization'
            write(luout, '(a)') '         effects can be removed by damping the electric field on polarizable sites'
            write(luout, '(a/)') '         (please see the manual for damping options).'
        end if
    end if

 10 format(12x,a16,5x,f20.12)
 20 format(/7x,a)
    write(luout, '(8x,a)') 'Polarizable embedding energy contributions:'
    write(luout, '(7x,a)') '---------------------------------------------'
    if (any(lmul)) then
        write(luout, 20) 'Electrostatic contributions:'
        write(luout, 10) 'Electronic      ', Ees(1,1)
        write(luout, 10) 'Nuclear         ', Ees(2,1)
        write(luout, 10) 'Multipole       ', Ees(3,1)
        write(luout, 10) 'Total           ', sum(Ees(:,1))
    end if
    if (any(lpol)) then
        write(luout, 20) 'Polarization contributions:'
        write(luout, 10) 'Electronic      ', Epol(1,1)
        write(luout, 10) 'Nuclear         ', Epol(2,1)
        if (any(lmul)) then
            write(luout, 10) 'Multipole       ', Epol(3,1)
        end if
        write(luout, 10) 'Total           ', sum(Epol(:,1))
    end if
    if (lvdw) then
        write(luout, 20) 'LJ energy contribution:'
        write(luout, 10) 'Dispersion 6    ', Elj(1)
        write(luout, 10) 'Repulsion 12    ', Elj(2)
    end if
    write(luout,'(/6x,a18,6x,f20.12)') 'Total PE energy: ', Epe(1)
    if (lpol(1,1) .and. pe_verbose) then
        write(luout, '(/15x,a)') 'Polarizable embedding information:'
        write(luout, '(14x,a/)') '------------------------------------'
        write(luout, '(23x,a)') 'Induced dipole moments'
        write(luout, '(7x,a,10x,a,14x,a,14x,a)') 'site', 'X', 'Y', 'Z'
        j = 1
        do i = 1, nsites
            if (zeroalphas(i)) cycle
            write(luout, '(5x,i6,3f15.8)') i, Mkinds(j:j+2,1)
            j = j + 3
        end do
        indtot = 0.0_dp
        j = 1
        do i = 1, npols
            indtot(1) = indtot(1) + Mkinds(j,1)
            indtot(2) = indtot(2) + Mkinds(j+1,1)
            indtot(3) = indtot(3) + Mkinds(j+2,1)
            j = j + 3
        end do
        write(luout, '(/22x,a)') 'Total induced dipole moment'
        write(luout, '(21x,a,14x,a,14x,a)') 'X', 'Y', 'Z'
        write(luout, '(11x,3f15.8/)') indtot
    end if

    if (lpol(1,1)) deallocate(Mkinds)

end subroutine print_info

!-----------------------------------------------------------------------------

subroutine electrostatic_operator(denmats, fckmats)

    use pe_variables

    real(dp), dimension(:), intent(in) :: denmats
    real(dp), dimension(:), intent(inout), optional :: fckmats

    integer :: i

    if (fock) then
        call es_multipoles(denmats, fckmats)
    else if (energy) then
        call es_multipoles(denmats)
    end if

end subroutine electrostatic_operator

!------------------------------------------------------------------------------

subroutine es_multipoles(denmats, fckmats)

    use pe_mpi
    use pe_variables
    use pe_utils, only: openfile
    use pe_blas_interfaces, only: dot
    use pe_multipole_tools, only: prefactors, Tk_tensor

    real(dp), dimension(:), intent(in) :: denmats
    real(dp), dimension(:), intent(inout), optional :: fckmats

    logical :: lexist
    integer :: lu
    integer :: site, ncomps
    integer :: i, j, k
    real(dp), dimension(3) :: Rsm, Rss
    real(dp), dimension(:), allocatable :: Tsm, factors
    real(dp), dimension(:,:), allocatable :: Mk_ints
    real(dp), dimension(:), allocatable :: tmpfckmats

    Ees = 0.0_dp

    if (myid == master) then
        inquire(file='pe_electrostatics.bin', exist=lexist)
    end if
#if defined(VAR_MPI)
    if (nprocs > 1) then
        call mpi_bcast(lexist, 1, lmpi, master, comm, ierr)
    end if
#endif
    if (lexist .and. fock .and. ((scf_cycle > 1) .or. pe_restart)) then
        if (myid == master) then
            call openfile('pe_electrostatics.bin', lu, 'old', 'unformatted')
            rewind(lu)
            read(lu) Ees(2,1), Ees(3,1), fckmats
            close(lu)
            do i = 1, ndens
                j = (i - 1) * nnbas + 1
                k = i * nnbas
                Ees(1,i) = dot(denmats(j:k), fckmats(j:k))
            end do
        end if
    else
        allocate(Mk_ints(nnbas,1))
        do site = site_start, site_finish
            ! electron - multipole interaction energy
            Mk_ints = 0.0_dp
            if (lmul(0)) then
                if (maxval(abs(M0s(:,site))) >= zero) then
                    call Mk_integrals('potential_derivative', Mk_ints, Rs(:,site), M0s(:,site))
                end if
            end if
            if (lmul(1)) then
                if (maxval(abs(M1s(:,site))) >= zero) then
                    call Mk_integrals('potential_derivative', Mk_ints, Rs(:,site), M1s(:,site))
                end if
            end if
            if (lmul(2)) then
                if (maxval(abs(M2s(:,site))) >= zero) then
                    call Mk_integrals('potential_derivative', Mk_ints, Rs(:,site), M2s(:,site))
                end if
            end if
            if (lmul(3)) then
                if (maxval(abs(M3s(:,site))) >= zero) then
                    call Mk_integrals('potential_derivative', Mk_ints, Rs(:,site), M3s(:,site))
                end if
            end if
            if (lmul(4)) then
                if (maxval(abs(M4s(:,site))) >= zero) then
                    call Mk_integrals('potential_derivative', Mk_ints, Rs(:,site), M4s(:,site))
                end if
            end if
            if (lmul(5)) then
                if (maxval(abs(M5s(:,site))) >= zero) then
                    call Mk_integrals('potential_derivative', Mk_ints, Rs(:,site), M5s(:,site))
                end if
            end if
            do i = 1, ndens
                j = (i - 1) * nnbas + 1
                k = i * nnbas
                Ees(1,i) = Ees(1,i) + dot(denmats(j:k), Mk_ints(:,1))
                if (fock) fckmats(j:k) = fckmats(j:k) + Mk_ints(:,1)
            end do

            ! nuclei - multipole interaction energy
            if (lmul(0)) then
                ncomps = size(M0s(:,site), 1)
                allocate(Tsm(ncomps))
                allocate(factors(ncomps))
                do i = 1, nnucs
                    Rsm = Rm(:,i) - Rs(:,site)
                    call Tk_tensor(Tsm, Rsm)
                    call prefactors(factors)
                    do j = 1, ncomps
                        Ees(2,1) = Ees(2,1) + factors(j) * M0s(j,site) * Tsm(j) * Zm(1,i)
                    end do
                end do
                deallocate(Tsm)
                deallocate(factors)
            end if
            if (lmul(1)) then
                ncomps = size(M1s(:,site), 1)
                allocate(Tsm(ncomps))
                allocate(factors(ncomps))
                do i = 1, nnucs
                    Rsm = Rm(:,i) - Rs(:,site)
                    call Tk_tensor(Tsm, Rsm)
                    call prefactors(factors)
                    do j = 1, ncomps
                        Ees(2,1) = Ees(2,1) + factors(j) * M1s(j,site) * Tsm(j) * Zm(1,i)
                    end do
                end do
                deallocate(Tsm)
                deallocate(factors)
            end if
            if (lmul(2)) then
                ncomps = size(M2s(:,site), 1)
                allocate(Tsm(ncomps))
                allocate(factors(ncomps))
                do i = 1, nnucs
                    Rsm = Rm(:,i) - Rs(:,site)
                    call Tk_tensor(Tsm, Rsm)
                    call prefactors(factors)
                    do j = 1, ncomps
                        Ees(2,1) = Ees(2,1) + factors(j) * M2s(j,site) * Tsm(j) * Zm(1,i)
                    end do
                end do
                deallocate(Tsm)
                deallocate(factors)
            end if
            if (lmul(3)) then
                ncomps = size(M3s(:,site), 1)
                allocate(Tsm(ncomps))
                allocate(factors(ncomps))
                do i = 1, nnucs
                    Rsm = Rm(:,i) - Rs(:,site)
                    call Tk_tensor(Tsm, Rsm)
                    call prefactors(factors)
                    do j = 1, ncomps
                        Ees(2,1) = Ees(2,1) + factors(j) * M3s(j,site) * Tsm(j) * Zm(1,i)
                    end do
                end do
                deallocate(Tsm)
                deallocate(factors)
            end if
            if (lmul(4)) then
                ncomps = size(M4s(:,site), 1)
                allocate(Tsm(ncomps))
                allocate(factors(ncomps))
                do i = 1, nnucs
                    Rsm = Rm(:,i) - Rs(:,site)
                    call Tk_tensor(Tsm, Rsm)
                    call prefactors(factors)
                    do j = 1, ncomps
                        Ees(2,1) = Ees(2,1) + factors(j) * M4s(j,site) * Tsm(j) * Zm(1,i)
                    end do
                end do
                deallocate(Tsm)
                deallocate(factors)
            end if
            if (lmul(5)) then
                ncomps = size(M5s(:,site), 1)
                allocate(Tsm(ncomps))
                allocate(factors(ncomps))
                do i = 1, nnucs
                    Rsm = Rm(:,i) - Rs(:,site)
                    call Tk_tensor(Tsm, Rsm)
                    call prefactors(factors)
                    do j = 1, ncomps
                        Ees(2,1) = Ees(2,1) + factors(j) * M5s(j,site) * Tsm(j) * Zm(1,i)
                    end do
                end do
                deallocate(Tsm)
                deallocate(factors)
            end if

            ! multipole - multipole interaction energy
            !do i = 1, nsites
            !    Rss = Rs(:,i) - Rs(:,site)
            !    if (lmul(0)) then
            !        ncomps = size(M0s(:,site), 1)
            !        allocate(Tsm(ncomps))
            !        allocate(factors(ncomps))
            !        call Tk_tensor(Tsm, Rsm)
            !        call prefactors(factors)
            !        do j = 1, ncomps
            !            Ees(3,1) = Ees(3,1) + factors(j) * M0s(j,site) * Tsm(j) * M0s(j,i)
            !        end do
            !        deallocate(Tsm)
            !        deallocate(factors)
            !    end if
            !end do
        end do
        deallocate(Mk_ints)

#if defined(VAR_MPI)
        if (myid == master .and. nprocs > 1) then
            call mpi_reduce(mpi_in_place, Ees(1,1), 3 * ndens, rmpi, mpi_sum, master, comm, ierr)
        else if (myid /= master) then
            call mpi_reduce(Ees(1,1), dummy_real, 3 * ndens, rmpi, mpi_sum, master, comm, ierr)
        end if
#endif

        if (fock) then
#if defined(VAR_MPI)
            if (myid == master .and. nprocs > 1) then
                allocate(tmpfckmats(ndens * nnbas))
                tmpfckmats = fckmats
                call mpi_reduce(mpi_in_place, fckmats(1), ndens * nnbas, rmpi, mpi_sum, master, comm, ierr)
            else if (myid /= master) then
                call mpi_reduce(fckmats(1), dummy_real, ndens * nnbas, rmpi, mpi_sum, master, comm, ierr)
            end if
#endif
            if (myid == master) then
                call openfile('pe_electrostatics.bin', lu, 'unknown', 'unformatted')
                rewind(lu)
                write(lu) Ees(2,1), Ees(3,1), fckmats
                close(lu)
            end if
#if defined(VAR_MPI)
            if (myid == master .and. nprocs > 1) then
                fckmats = tmpfckmats
                deallocate(tmpfckmats)
            end if
#endif
        end if
    end if

end subroutine es_multipoles

!------------------------------------------------------------------------------

subroutine polarization_operator(denmats, fckmats)

    use pe_mpi
    use pe_variables
    use pe_potential_derivatives
    use pe_integral_interfaces
    use pe_blas_interfaces, only: dot
    use pe_induced_moments, only: induced_moments

    real(dp), dimension(:), intent(in), optional :: denmats
    real(dp), dimension(:), intent(inout), optional :: fckmats

    integer :: site
    integer :: i, j, k, l, m
    real(dp), dimension(3) :: indtot
    real(dp), dimension(:,:), allocatable :: Vels, Vtots, Qinds
    real(dp), dimension(:), allocatable :: Vnucs, Vmuls, Vfds
    real(dp), dimension(:,:,:), allocatable :: Vel_ints
    real(dp), dimension(:,:), allocatable :: Fels, Ftots, Mkinds
    real(dp), dimension(:), allocatable :: Fnucs, Fmuls, Fext
    real(dp), dimension(:,:,:), allocatable :: Fel_ints

    if (lpol(1,1)) then
        allocate(Mkinds(3*npols,ndens), Fels(3*npols,ndens), Ftots(3*npols,ndens))
        Fels = 0.0_dp; Ftots = 0.0_dp
        allocate(Fnucs(3*npols), Fmuls(3*npols), Fext(3*npols))
        Fnucs = 0.0_dp; Fmuls = 0.0_dp; Fext= 0.0_dp
    end if

    if (response) then
        if (lpol(1,1)) then
            call electron_fields(Fels, denmats)
            if (myid == master) then
                do i = 1, ndens
                    Ftots(:,i) = Fels(:,i)
                end do
            end if
        end if
    else
        if (lpol(1,1)) then
            if (pe_field) then
                do site = site_start, site_finish
                    if (zeroalphas(site)) cycle
                    do i = 1, 3
                        Fext(i+3*(site-1)) = efields(i)
                    end do
                end do
            end if
            call electron_fields(Fels, denmats)
            call nuclear_fields(Fnucs)
            call multipole_fields(Fmuls)
            if (myid == master) then
                do i = 1, ndens
                    Ftots(:,i) = Fels(:,i) + Fnucs + Fmuls + Fext
                end do
            end if
        end if
    end if
    if (any(lpol)) then
        call induced_moments(Mkinds=Mkinds, Fks=Ftots)
    end if
    if (.not. (response) .and. myid == master) then
        do i = 1, ndens
            if (any(lpol)) then
                if (lpol(1,1)) then
                    Epol(1,i) = - 0.5_dp * dot(Mkinds(:,i), Fels(:,i))
                    Epol(2,i) = - 0.5_dp * dot(Mkinds(:,i), Fnucs)
                    if (any(lmul)) Epol(3,i) = - 0.5_dp * dot(Mkinds(:,i), Fmuls)
                end if
            end if
        end do
    end if
    if (myid == master .and. (pe_debug .or. pe_verbose) .and. .not. energy) then
        write(luout, '(/6x,a)') 'Info from polarizable embedding library:'
        write(luout, '(5x,a/)') '------------------------------------------'
        do i = 1, ndens
            write(luout, '(7x,a,i3)') 'Input density no.: ', i
            if (lpol(1,1) .and. .not. pe_debug) then
                indtot = 0.0_dp
                k = 1
                do j = 1, npols
                    indtot(1) = indtot(1) + Mkinds(k,i)
                    indtot(2) = indtot(2) + Mkinds(k+1,i)
                    indtot(3) = indtot(3) + Mkinds(k+2,i)
                    k = k + 3
                end do
                write(luout, '(/24x,a)') 'Total induced dipole moment'
                write(luout, '(23x,a,14x,a,14x,a)') 'X', 'Y', 'Z'
                write(luout, '(13x,3f15.8/)') indtot
            end if
            if (any(lpol) .and. pe_debug) then
                write(luout, '(/15x,a)') 'Total electric field at polarizable sites:'
                write(luout, '(9x,a,10x,a,14x,a,14x,a)') 'site', 'X', 'Y', 'Z'
                k = 1
                do site = 1, nsites
                    if (zeroalphas(site)) cycle
                    write(luout, '(7x,i6,3f15.8)') site, Ftots(k:k+2,i)
                    k = k + 3
                end do
                write(luout, '(/25x,a)') 'Induced dipole moments'
                write(luout, '(9x,a,10x,a,14x,a,14x,a)') 'site', 'X', 'Y', 'Z'
                k = 1
                do site = 1, nsites
                    if (zeroalphas(site)) cycle
                    write(luout, '(7x,i6,3f15.8)') site, Mkinds(k:k+2,i)
                    k = k + 3
                end do
            end if
        end do
    end if
#if defined(VAR_MPI)
    if (myid == master .and. nprocs > 1) then
        do i = 1, ndens
            if (lpol(1,1)) then
                displs(0) = 0
                do j = 1, nprocs
                    displs(j) = displs(j-1) + poldists(j-1)
                end do
                call mpi_scatterv(Mkinds(1,i), poldists, displs, rmpi, mpi_in_place, dummy_int, rmpi, master, comm, ierr)
            end if
        end do
    else if (myid /= master) then
        do i = 1, ndens
            if (lpol(1,1)) then
                call mpi_scatterv(dummy_real, poldists, displs, rmpi, Mkinds(1,i), poldists(myid), rmpi, master, comm, ierr)
            end if
        end do
    end if
#endif
    if (fock .or. response) then
        if (any(lpol)) then
            allocate(Fel_ints(nnbas,3,1))
            i = 0
            do site = site_start, site_finish
                if (zeroalphas(site)) cycle
                call Tk_integrals('potential_derivative', Rs(:,site), Fel_ints)
                do j = 1, 3
                    do k = 1, ndens
                        l = (k - 1) * nnbas + 1
                        m = k * nnbas
                        fckmats(l:m) = fckmats(l:m) - Mkinds(i+j,k) * Fel_ints(:,j,1)
                    end do
                end do
                i = i + 3
            end do
            deallocate(Fel_ints)
        end if
    end if

    if (any(lpol)) then
        deallocate(Mkinds, Fels, Ftots)
        deallocate(Fnucs, Fmuls)
    end if

end subroutine polarization_operator

!------------------------------------------------------------------------------
!> Computes the effective external field contributions to the electric dipole
! operator in the polarizable embedding scheme.
!! @author Nanna Holmgaard List

subroutine get_effdipole(fckmats)

    use pe_mpi
    use pe_variables
    use pe_integral_interfaces
    use pe_blas_interfaces, only: dot
    use pe_induced_moments, only: induced_moments
    use pe_multipole_tools, only: multipole_derivative

    real(dp), dimension(3) :: Rsc, Fp
    real(dp), dimension(3,3,nnucs+1+ncrds) :: lf_tensor
    real(dp), dimension(3,3) :: alpha_mol
    real(dp), dimension(:), intent(out) :: fckmats
    real(dp), dimension(:,:), allocatable :: Fext, Mkinds, Binv_up
    real(dp), dimension(:,:), allocatable :: Minds_int
    real(dp), dimension(:,:,:), allocatable :: Fel_ints, Fp_tot
    real(dp), dimension(:), allocatable :: B, Binv
    logical, parameter :: local_debug = .false.
    integer :: i, l, m, j, k, ioff, site, info

    if (.not. lpol(1,1)) then
        stop 'ERROR: non-polarizable environment, i.e. no field factors'
    end if

    allocate(Fext(3*npols,3),Mkinds(3*npols,3),Fel_ints(nnbas,3,1))

    Fext = 0.0_dp; Fel_ints = 0.0_dp; Mkinds = 0.0_dp
    l = 1
    do site = 1, nsites
        if (zeroalphas(site)) then
           cycle
        end if
        Fext(l,1) = 1.0_dp
        Fext(l+1,2) = 1.0_dp
        Fext(l+2,3) = 1.0_dp
        l = l + 3
    end do

    call induced_moments(Mkinds, Fext)

    if (local_debug) then
        write(luout, '(/15x,a)') 'Local field factor information:'
        write(luout, '(14x,a/)') '------------------------------------'
        do m = 1, 3
            j = 1
            write(luout, '(23x,a,i1)') 'Derivatives of induced dipole moments wrt. Fext_m: ', m
            write(luout, '(7x,a,10x,a,14x,a,14x,a)') 'site', 'X', 'Y', 'Z'
            do i = 1, nsites
                if (zeroalphas(i)) cycle
                write(luout, '(5x,i6,3f15.8)') i, Mkinds(j:j+2,m)
                j = j + 3
            end do
        end do
    end if

#if defined(VAR_MPI)
    if (myid == master .and. nprocs > 1) then
        do j = 1, 3
            displs(0) = 0
            do i = 1, nprocs
                displs(i) = displs(i-1) + poldists(i-1)
            end do
            call mpi_scatterv(Mkinds(1,j), poldists, displs, rmpi, mpi_in_place, dummy_int, rmpi, master, comm, ierr)
        end do
    else if (myid /= master) then
        do j = 1, 3
            call mpi_scatterv(dummy_real, poldists, displs, rmpi, Mkinds(1,j), poldists(myid), rmpi, master, comm, ierr)
        end do
    end if
#endif

    l = 0
    do site = site_start, site_finish
        if (zeroalphas(site)) cycle
        call Tk_integrals('potential_derivative', Rs(:,site), Fel_ints)
        do m = 1, 3 !x, y, z dipole components
            do i = 1, 3 !x, y, z external field
                j = (i - 1) * nnbas + 1
                k = i * nnbas
                fckmats(j:k) = fckmats(j:k) - Mkinds(l+m,i) * Fel_ints(:,m,1)
            end do
        end do
        l = l + 3
    end do

#if defined(VAR_MPI)
    if (myid == master .and. nprocs > 1) then
        call mpi_reduce(mpi_in_place, fckmats(1), 3 * nnbas , rmpi, mpi_sum, master, comm, ierr)
    else if (myid /= master) then
        call mpi_reduce(fckmats(1), dummy_real, 3 * nnbas, rmpi, mpi_sum, master, comm, ierr)
    end if
#endif

    lf_tensor = 0.0_dp
    alpha_mol = 0.0_dp
    l = 1
    do site = site_start, site_finish
        if (zeroalphas(site)) cycle
!       lf_tensor at atomic sites
        do j = 1, nnucs
           Rsc = Rm(:,j) - Rs(:,site)
           do m = 1, 3
              Fp = 0.0_dp
              call multipole_derivative(Fp, Rsc, Mkinds(l:l+2,m))
              lf_tensor(:,m,j) = lf_tensor(:,m,j) + Fp
           end do
        end do
!       lf_tensor at COM
        Rsc = core_com(:) - Rs(:,site)
        do m = 1, 3
            Fp = 0.0_dp
            call multipole_derivative(Fp, Rsc, Mkinds(l:l+2,m))
            lf_tensor(:,m,nnucs+1) = lf_tensor(:,m,nnucs+1) + Fp
        end do
!       lf_tensor at external sites given in input
        if (allocated(crds)) then
            do j = 1, ncrds
               Rsc = crds(:,j) - Rs(:,site)
               do m = 1, 3
                  Fp = 0.0_dp
                  call multipole_derivative(Fp, Rsc, Mkinds(l:l+2,m))
                  lf_tensor(:,m,nnucs+1+j) = lf_tensor(:,m,nnucs+1+j) + Fp
               end do
            end do
        end if
        l = l + 3
    end do

#if defined(VAR_MPI)
    if (myid == master .and. nprocs > 1) then
        do m = 1, 3
            call mpi_reduce(mpi_in_place, alpha_mol(1,m), 3, rmpi, mpi_sum, master, comm, ierr)
            do j = 1, nnucs + 1 + ncrds
                call mpi_reduce(mpi_in_place, lf_tensor(1,m,j), 3, rmpi, mpi_sum, master, comm, ierr)
            end do
        end do
    else if (myid /= master) then
        do m = 1, 3
            call mpi_reduce(alpha_mol(1,m), dummy_real, 3, rmpi, mpi_sum, master, comm, ierr)
            do j = 1, nnucs + 1 + ncrds
                call mpi_reduce(lf_tensor(1,m,j), dummy_real, 3, rmpi, mpi_sum, master, comm, ierr)
            end do
        end do
    end if
#endif

    if (myid == master) then
        do j = 1, nnucs + 1 + ncrds
           do m = 1, 3
               lf_tensor(m,m,j) = lf_tensor(m,m,j) + 1.0_dp
           end do
        end do
    end if

    if (myid == master) then
        write(luout, '(/9x,a)') 'Effective external field (EEF) information:'
        write(luout, '(9x,a/)') '-------------------------------------------'

        do j = 1, nnucs + 1 + ncrds
           if (j > nnucs + 1) then
               write(luout, '(20x,a,i5)') 'EEF tensor at site: ', j
           else if (j == nnucs + 1) then
               write(luout, '(21x,a)') 'EEF tensor at COM'
           else
               cycle
!               write(luout, '(/20x,a,5i)') 'EEF tensor at nuclei: ', j
           end if
           do m = 1, 3
               write(luout, '(5x,3f15.8)') (lf_tensor(i,m,j), i = 1,3)
           end do
        end do
        do j = 1, 3
            do m = 1, 3
                alpha_mol(j,m) = dot(Fext(:,j),Mkinds(:,m))
            end do
        end do
        write(luout, '(/13x,a)') 'Environment polarizability (a.u.)'
        do m = 1, 3
            write(luout, '(5x,3f15.8)') (alpha_mol(m,i), i = 1,3)
        end do
    end if 
    deallocate(Fext)
    deallocate(Mkinds)
    deallocate(Fel_ints)

end subroutine get_effdipole

!------------------------------------------------------------------------------

subroutine get_molecular_gradient(denmats, fckmats, molgrads)

    use pe_mpi
    use pe_variables

    real(dp), dimension(:), intent(in), optional :: denmats
    real(dp), dimension(:), intent(out), optional :: fckmats
    real(dp), dimension(:), intent(out), optional :: molgrads

    if (present(molgrads) .and. .not. present(denmats)) then
        stop 'ERROR: missing denmats in pe_molecular_gradients'
    else if (present(denmats) .and. .not. present(molgrads)) then
        stop 'ERROR: missing molgrads in pe_molecular_gradients'
    else if (.not. present(molgrads) .and. .not. present(fckmats)) then
        stop 'ERROR: nothing to do in pe_molecular_gradients'
    end if

    if (any(lmul)) then
        if (present(molgrads) .and. present(fckmats)) then
            call es_gradients(denmats=denmats, fckmats=fckmats, molgrads=molgrads)
        else if (present(molgrads)) then
            call es_gradients(denmats=denmats, molgrads=molgrads)
        else if (present(fckmats)) then
            call es_gradients(fckmats=fckmats)
        end if
    end if
    ! true when having (iso)alphas
    if (lpol(1,1)) then
        if (present(molgrads) .and. present(fckmats)) then
            call pol_gradients(denmats=denmats, fckmats=fckmats, molgrads=molgrads)
        else if (present(molgrads)) then
            call pol_gradients(denmats=denmats, molgrads=molgrads)
        else if (present(fckmats)) then
            call pol_gradients(fckmats=fckmats)
        end if
    end if
    if (present(molgrads)) then
        if (lvdw) then
            call lj_gradients(molgrads)
        else
           write(luout, *) 'WARNING: No LJ contribution to geometrical derivatives is included'
        end if
    end if

#if defined(VAR_MPI)
    if (myid == master .and. nprocs > 1) then
        call mpi_reduce(mpi_in_place, molgrads(1), 3 * nnucs , rmpi, mpi_sum, master, comm, ierr)
    else if (myid /= master) then
        call mpi_reduce(molgrads(1), dummy_real, 3 * nnucs, rmpi, mpi_sum, master, comm, ierr)
    end if
#endif

end subroutine get_molecular_gradient

!------------------------------------------------------------------------------

subroutine pol_gradients(denmats, fckmats, molgrads)

    use pe_mpi
    use pe_variables
    use pe_utils, only: openfile
    use pe_blas_interfaces, only: dot, spmv
    use pe_multipole_tools, only: multipole_derivative
    use pe_potential_derivatives
    use pe_induced_moments, only: induced_moments

    real(dp), dimension(:), intent(in), optional :: denmats
    real(dp), dimension(:), intent(inout), optional :: molgrads
    real(dp), dimension(:), intent(inout), optional :: fckmats

    real(dp), dimension(3) :: Rms
    real(dp), dimension(6) :: grdFnucs
    real(dp), dimension(:), allocatable :: Fnucs, Fmuls
    real(dp), dimension(:,:), allocatable :: Mkinds, Fktots, Fels
    real(dp), dimension(:,:), allocatable :: Mk_ints

    integer :: i, j, k, l, site, lu
    logical :: lexist

    if (present(molgrads) .and. .not. present(denmats)) then
        stop 'ERROR: missing denmats in pol_gradients'
    else if (present(denmats) .and. .not. present(molgrads)) then
        stop 'ERROR: missing molgrads in pol_gradients'
    else if (.not. present(molgrads) .and. .not. present(fckmats)) then
        stop 'ERROR: nothing to do in pol_gradients'
    end if

    allocate(Mk_ints(nnbas,3*nnucs))
    allocate(Mkinds(3*npols,ndens))

    Mkinds = 0.0_dp

!   Read in induced moments
    if (myid == master) then
        inquire(file='pe_induced_moments.bin', exist=lexist)
        if (.not. lexist) then
            stop 'ERROR: pe_induced_moments.bin does not exist in pol_gradients'
        end if
    end if

#if defined(VAR_MPI)
    if (nprocs > 1) then
        call mpi_bcast(lexist, 1, lmpi, master, comm, ierr)
    end if
#endif

    if (lexist) then
        if (myid == master) then
            call openfile('pe_induced_moments.bin', lu, 'old', 'unformatted')
            rewind(lu)
            read(lu) Mkinds
            close(lu)
        end if
    else if (present(fckmats)) then
        stop 'ERROR: Updating induced dipoles according to a wrong density matrix'
    else
        allocate(Fels(3*npols,ndens), Fktots(3*npols,ndens))
        allocate(Fnucs(3*npols), Fmuls(3*npols))
        Fels = 0.0_dp; Fktots = 0.0_dp
        Fnucs = 0.0_dp; Fmuls = 0.0_dp
        call electron_fields(Fels, denmats)
        call nuclear_fields(Fnucs)
        call multipole_fields(Fmuls)
        do i = 1, ndens
            Fktots(1:3*npols,i) = Fels(:,i) + Fnucs + Fmuls
        end do
        call induced_moments(Mkinds, Fktots)
        deallocate(Fels, Fktots, Fnucs, Fmuls)
    end if

#if defined(VAR_MPI)
    if (myid == master .and. nprocs > 1) then
        displs(0) = 0
        do i = 1, nprocs
            displs(i) = displs(i-1) + poldists(i-1)
        end do
        call mpi_scatterv(Mkinds(1,1), poldists, displs, rmpi, mpi_in_place, dummy_int, rmpi, master, comm, ierr)
    else if (myid /= master) then
        call mpi_scatterv(dummy_real, poldists, displs, rmpi, Mkinds(1,1), poldists(myid), rmpi, master, comm, ierr)
    end if
#endif

    ! 1st derivatives of nuclear and electron fields multiplied with induced dipoles
    l = 1
    do site = site_start, site_finish
        if (zeroalphas(site)) cycle
        if (present(molgrads)) then
            k = 1
            do i = 1, nnucs
                grdFnucs = 0.0_dp
                Rms = Rs(:,site) - Rm(:,i)
                call multipole_derivative(grdFnucs, Rms, Zm(:,i))
                call spmv(grdFnucs, Mkinds(l:l+2,1), molgrads(k:k+2), 'L', -1.0_dp, 1.0_dp)
                k = k + 3
            end do
        end if

        Mk_ints = 0.0_dp
        call Mk_integrals('molecular_gradient', Mk_ints, Rs(:,site), Mkinds(l:l+2,1))
        l = l + 3

        if (present(molgrads)) then
            do i = 1, 3 * nnucs
                molgrads(i) = molgrads(i) + dot(denmats(1:nnbas), Mk_ints(:,i))
            end do
        end if
        if (present(fckmats)) then
            do i = 1, 3 * nnucs
                j = (i - 1) * nnbas + 1
                k = i * nnbas
                fckmats(j:k) = fckmats(j:k) + Mk_ints(:,i)
            end do
        end if
    end do

    deallocate(Mkinds, Mk_ints)

end subroutine pol_gradients

!------------------------------------------------------------------------------

subroutine es_gradients(denmats, fckmats, molgrads)

!   Calculates gradient contribution from QM nuclear - multipole
!   and QM electron - multipole interaction energies

    use pe_variables
    use pe_blas_interfaces, only: dot
    use pe_multipole_tools, only: multipole_derivative

    real(dp), dimension(:), intent(in), optional :: denmats
    real(dp), dimension(:), intent(inout), optional :: fckmats
    real(dp), dimension(:), intent(inout), optional :: molgrads

    integer :: i, j, k, site
    real(dp), dimension(3) :: Rsm
    real(dp), dimension(3) :: mul_grad
    real(dp), dimension(:,:), allocatable :: Mk_ints

    if (present(molgrads) .and. .not. present(denmats)) then
        stop 'ERROR: missing denmats in es_gradients'
    else if (present(denmats) .and. .not. present(molgrads)) then
        stop 'ERROR: missing molgrads in es_gradients'
    else if (.not. present(molgrads) .and. .not. present(fckmats)) then
        stop 'ERROR: nothing to do in es_gradients'
    end if

    allocate(Mk_ints(nnbas,3*nnucs))

    Mk_ints = 0.0_dp
    do site = site_start, site_finish
        ! 1st derivative of nuclei - multipole interaction energy
        if (present(molgrads)) then
            k = 1
            do i = 1, nnucs
                mul_grad = 0.0_dp
                Rsm = Rm(:,i) - Rs(:,site)
                if (lmul(0)) then
                    if (maxval(abs(M0s(:,site))) >= zero) then
                        call multipole_derivative(mul_grad, Rsm, M0s(:,site))
                    end if
                end if
                if (lmul(1)) then
                    if (maxval(abs(M1s(:,site))) >= zero) then
                        call multipole_derivative(mul_grad, Rsm, M1s(:,site))
                    end if
                end if
                if (lmul(2)) then
                    if (maxval(abs(M2s(:,site))) >= zero) then
                        call multipole_derivative(mul_grad, Rsm, M2s(:,site))
                    end if
                end if
                if (lmul(3)) then
                    if (maxval(abs(M3s(:,site))) >= zero) then
                        call multipole_derivative(mul_grad, Rsm, M3s(:,site))
                    end if
                end if
                if (lmul(4)) then
                    if (maxval(abs(M4s(:,site))) >= zero) then
                        call multipole_derivative(mul_grad, Rsm, M4s(:,site))
                    end if
                end if
                if (lmul(5)) then
                    if (maxval(abs(M5s(:,site))) >= zero) then
                        call multipole_derivative(mul_grad, Rsm, M5s(:,site))
                    end if
                end if
                ! minus here because multipole_derivative gives opposite sign
                molgrads(k:k+2) = molgrads(k:k+2) - Zm(1,i) * mul_grad
                k = k + 3
            end do
        end if

        ! 1st derivative of electron - multipole interaction energy
        if (lmul(0)) then
            if (maxval(abs(M0s(:,site))) >= zero) then
                call Mk_integrals('molecular_gradient', Mk_ints, Rs(:,site), M0s(:,site))
            end if
        end if
        if (lmul(1)) then
            if (maxval(abs(M1s(:,site))) >= zero) then
                call Mk_integrals('molecular_gradient', Mk_ints, Rs(:,site), M1s(:,site))
            end if
        end if
        if (lmul(2)) then
            if (maxval(abs(M2s(:,site))) >= zero) then
                call Mk_integrals('molecular_gradient', Mk_ints, Rs(:,site), M2s(:,site))
            end if
        end if
        if (lmul(3)) then
            if (maxval(abs(M3s(:,site))) >= zero) then
                call Mk_integrals('molecular_gradient', Mk_ints, Rs(:,site), M3s(:,site))
            end if
        end if
        if (lmul(4)) then
            if (maxval(abs(M4s(:,site))) >= zero) then
                call Mk_integrals('molecular_gradient', Mk_ints, Rs(:,site), M4s(:,site))
            end if
        end if
        if (lmul(5)) then
            if (maxval(abs(M5s(:,site))) >= zero) then
                call Mk_integrals('molecular_gradient', Mk_ints, Rs(:,site), M5s(:,site))
            end if
        end if
    end do
    if (present(molgrads)) then
        do i = 1, 3 * nnucs
            molgrads(i) = molgrads(i) + dot(denmats(1:nnbas), Mk_ints(:,i))
        end do
    end if
    if (present(fckmats)) then
        do i = 1, 3 * nnucs
            j = (i - 1) * nnbas + 1
            k = i * nnbas
            fckmats(j:k) = fckmats(j:k) + Mk_ints(:,i)
        end do
    end if

    deallocate(Mk_ints)

end subroutine es_gradients

!------------------------------------------------------------------------------

subroutine lj_energy

    use pe_mpi
    use pe_variables
    use pe_constants
    use pe_blas_interfaces, only: nrm2

    integer :: i, site
    real(dp) :: eps_sm, r_sm, term_6, norm
    real(dp), dimension(3) :: Rsm

! TODO: Also be able to handle bond midpoints
    Elj    = 0.0_dp

    if (qmLJsites /= nnucs) stop 'Number of LJ sites does not match the&
                                 & number of nuclei in the core region.'

    do site = site_start, site_finish
        if (Zs(1,site) <= zero) cycle
        do i = 1, nnucs
            r_sm = 0.5_dp * (LJs(1,site) + qmLJs(1,i)) * aa2bohr
            eps_sm = sqrt(LJs(2,site) * qmLJs(2,i)) * kcal2hartree
            Rsm = Rm(:,i) - Rs(:,site)
            norm = nrm2(Rsm)
            term_6 = (r_sm / norm)**(6.0_dp)
            Elj(1) = Elj(1) - 2.0_dp * eps_sm * term_6
            Elj(2) = Elj(2) + eps_sm * term_6**(2.0_dp)
        end do
    end do

#if defined(VAR_MPI)
    if (myid == master .and. nprocs > 1) then
        call mpi_reduce(mpi_in_place, Elj(1), 2, rmpi, mpi_sum, master, comm, ierr)
    else if (myid /= master) then
        call mpi_reduce(Elj(1), dummy_real, 2, rmpi, mpi_sum, master, comm, ierr)
    end if
#endif

end subroutine lj_energy

!------------------------------------------------------------------------------

subroutine lj_gradients(molgrads)

    use pe_variables
    use pe_constants
    use pe_blas_interfaces, only: nrm2

    real(dp), dimension(:), intent(inout) :: molgrads

    integer :: i, j, site
    real(dp) :: eps_sm, r_sm, term_6, norm
    real(dp), dimension(3) :: fac, Rsm

    do site = site_start, site_finish
        if (Zs(1,site) <= zero) cycle
        j = 1
        do i = 1, nnucs
            r_sm = 0.5_dp * (LJs(1,site) + qmLJs(1,i)) * aa2bohr
            eps_sm = (sqrt(LJs(2,site) * qmLJs(2,i))) * kcal2hartree
            Rsm = Rm(:,i) - Rs(:,site)
            norm = nrm2(Rsm)
            fac = - 12.0_dp * eps_sm * Rsm / (norm**(2.0_dp))
            term_6 = (r_sm / norm)**(6.0_dp)
            molgrads(j:j+2) = molgrads(j:j+2) + fac * (term_6**(2.0_dp) - term_6)
            j = j + 3
        end do
    end do

end subroutine lj_gradients

!------------------------------------------------------------------------------

subroutine get_magnetic_gradient(fckmats)

    use pe_mpi
    use pe_variables
    use pe_utils, only: openfile

    real(dp), dimension(:), intent(out) :: fckmats

    integer :: i, lu
    real(dp), dimension(:,:), allocatable :: Mkinds
    logical :: lexist

    ! static multipole moment contribution
    if (lmul(0)) call lao_multipoles(M0s, fckmats)
    if (lmul(1)) call lao_multipoles(M1s, fckmats)
    if (lmul(2)) call lao_multipoles(M2s, fckmats)
    if (lmul(3)) call lao_multipoles(M3s, fckmats)
    if (lmul(4)) call lao_multipoles(M4s, fckmats)
    if (lmul(5)) call lao_multipoles(M5s, fckmats)

    ! polarization contribution
    if (lpol(1,1)) then
        allocate(Mkinds(3,npols))
        if (myid == master) then
            inquire(file='pe_induced_moments.bin', exist=lexist)
        end if
#if defined(VAR_MPI)
        if (nprocs > 1) then
            call mpi_bcast(lexist, 1, lmpi, master, comm, ierr)
        end if
#endif
        if (lexist) then
            if (myid == master) then
                call openfile('pe_induced_moments.bin', lu, 'old', 'unformatted')
                rewind(lu)
                read(lu) Mkinds
                close(lu)
            end if
        else
            stop 'ERROR: pe_induced_moments.bin does not exist'
        end if
#if defined(VAR_MPI)
        if (myid == master .and. nprocs > 1) then
            displs(0) = 0
            do i = 1, nprocs
                displs(i) = displs(i-1) + poldists(i-1)
            end do
            call mpi_scatterv(Mkinds(1,1), poldists, displs, rmpi, mpi_in_place, dummy_int, rmpi, master, comm, ierr)
        else if (myid /= master) then
            call mpi_scatterv(dummy_real, poldists, displs, rmpi, Mkinds(1,1), poldists(myid), rmpi, master, comm, ierr)
        end if
#endif
        call lao_induced_multipoles(Mkinds, fckmats)
        deallocate(Mkinds)
    endif

#if defined(VAR_MPI)
    if (myid == master .and. nprocs > 1) then
        call mpi_reduce(mpi_in_place, fckmats(1), 3 * nnbas, rmpi, mpi_sum, master, comm, ierr)
    else if (myid /= master) then
        call mpi_reduce(fckmats(1), dummy_real, 3 * nnbas, rmpi, mpi_sum, master, comm, ierr)
    end if
#endif

end subroutine get_magnetic_gradient

!------------------------------------------------------------------------------

subroutine lao_multipoles(Mks, fckmats)

    use pe_variables

    real(dp), dimension(:,:), intent(in) :: Mks
    real(dp), dimension(:), intent(inout) :: fckmats

    integer :: i, j, k
    integer :: site
    real(dp), dimension(:,:), allocatable :: Mk_ints

    allocate(Mk_ints(nnbas,3))

    Mk_ints = 0.0_dp
    do site = site_start, site_finish
        if (maxval(abs(Mks(:,site))) < zero) cycle
        call Mk_integrals('magnetic_gradient', Mk_ints, Rs(:,site), Mks(:,site))
    end do

    ! update for all directions of the magnetic field
    do i = 1, 3
        j = (i - 1) * nnbas + 1
        k = i * nnbas
        fckmats(j:k) = fckmats(j:k) + Mk_ints(:,i)
    end do

    deallocate(Mk_ints)

end subroutine lao_multipoles

!------------------------------------------------------------------------------

subroutine lao_induced_multipoles(Mks, fckmats)

    use pe_variables

    real(dp), dimension(:,:), intent(in) :: Mks
    real(dp), dimension(:), intent(inout) :: fckmats

    integer :: i, j, k
    integer :: site
    real(dp), dimension(:,:), allocatable :: Mk_ints

    allocate(Mk_ints(nnbas,3))

    Mk_ints = 0.0_dp
    i = 1
    do site = site_start, site_finish
        if (zeroalphas(site)) cycle
        call Mk_integrals('magnetic_gradient', Mk_ints, Rs(:,site), Mks(:,i))
        i = i + 1
    end do

    ! update for all directions of the magnetic field
    do i = 1, 3
        j = (i - 1) * nnbas + 1
        k = i * nnbas
        fckmats(j:k) = fckmats(j:k) + Mk_ints(:,i)
    end do

    deallocate(Mk_ints)

end subroutine lao_induced_multipoles

!------------------------------------------------------------------------------

subroutine Mk_integrals(inttype, Mk_ints, Rij, Mk)

    use pe_variables
    use pe_integral_interfaces
    use pe_multipole_tools, only: prefactors

    character(*), intent(in) :: inttype
    real(dp), dimension(:,:), intent(inout) :: Mk_ints
    real(dp), dimension(3), intent(in) :: Rij
    real(dp), dimension(:), intent(in) :: Mk

    integer :: i, j
    integer :: ncomps, nprops
    real(dp), dimension(:), allocatable :: factors
    real(dp), dimension(:,:,:), allocatable :: Tk_ints

    ncomps = size(Mk)

    if (inttype == 'potential_derivative') then
        nprops = 1
    else if (inttype == 'magnetic_gradient') then
        nprops = 3
    else if (inttype == 'molecular_gradient') then
        nprops = 3 * nnucs
    end if

    allocate(Tk_ints(nnbas,ncomps,nprops))

    call Tk_integrals(inttype, Rij, Tk_ints)

    allocate(factors(ncomps))
    call prefactors(factors)

    ! multiply T^(k) integrals with multipole to get M^(k) integrals
    do j = 1, nprops
        do i = 1, ncomps
            Mk_ints(:,j) = Mk_ints(:,j) + factors(i) * Mk(i) * Tk_ints(:,i,j)
        end do
    end do

    deallocate(factors)
    deallocate(Tk_ints)

end subroutine Mk_integrals

!------------------------------------------------------------------------------

end module pe_operators
