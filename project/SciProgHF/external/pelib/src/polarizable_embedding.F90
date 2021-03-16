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
module polarizable_embedding

    use pe_precision

    implicit none

    private

    ! public subroutines/functions
    public :: pe_init, pe_finalize, pe_input_reader, pe_master
#if defined(VAR_MPI)
    public :: pe_slave
#endif

contains

!------------------------------------------------------------------------------

subroutine pe_init(lupri, coords, charges)

    ! Initialization routine for the PE library.

    use pe_variables
    use pe_constants
    use pe_integral_interfaces
    use pe_blas_interfaces, only: nrm2
    use pe_utils, only: charge2mass, charge2elem

    integer :: lupri
    real(dp), dimension(:), intent(in), optional :: charges
    real(dp), dimension(:,:), intent(in), optional :: coords

    integer :: i, j, k, l
    integer :: idx, jdx, kdx, nidx
    integer, dimension(:), allocatable :: idxs
    logical, dimension(:), allocatable :: redists
    logical :: lexist, skip
    real(dp) :: rclose, mass, totmas, redist_charge
    real(dp), dimension(3) :: Rsm, Rji
    character(len=80) :: fd

    if (allocated(Rm) .and. allocated(Zm)) then
        Rm(:,:) = coords
        synced = .false.
        scf_cycle = 0
        core_com = 0.0_dp
        totmas = 0.0_dp
        do i = 1, nnucs
            mass = charge2mass(Zm(1,i))
            totmas = totmas + mass
            core_com = core_com + mass * Rm(:,i)
        end do
        core_com = core_com / totmas
        call init_pe_integral_interfaces(core_com, Rm)
        return
    end if

    luout = lupri

    initialized = .true.

    if (present(coords) .and. present(charges)) then
        nnucs = size(charges)
        allocate(Rm(3,nnucs), Zm(1,nnucs), core_com(3))
        Rm(:,:) = coords
        Zm(1,:) = charges
        core_com = 0.0_dp
        totmas = 0.0_dp
        do i = 1, nnucs
            mass = charge2mass(Zm(1,i))
            totmas = totmas + mass
            core_com = core_com + mass * Rm(:,i)
        end do
        core_com = core_com / totmas
        call init_pe_integral_interfaces(core_com, Rm)
    else if (present(coords) .and. .not. present(charges)) then
        stop 'ERROR in pe_init: coords present but charges missing'
    else if (.not. present(coords) .and. present(charges)) then
        stop 'ERROR in pe_init: charges present but coords missing'
    end if

    ! setting up grid for MEP and CUBE calculation
    if (pe_cube) then
        origin(1) = minval(Rm(1,:)) - xsize
        origin(2) = minval(Rm(2,:)) - ysize
        origin(3) = minval(Rm(3,:)) - zsize
        step(1) = 1.0_dp / xgrid
        step(2) = 1.0_dp / ygrid
        step(3) = 1.0_dp / zgrid
        xsteps = int((maxval(Rm(1,:)) + xsize - origin(1)) / step(1)) + 1
        ysteps = int((maxval(Rm(2,:)) + ysize - origin(2)) / step(2)) + 1
        zsteps = int((maxval(Rm(3,:)) + zsize - origin(3)) / step(3)) + 1
        npoints = 0
        npoints = xsteps * ysteps * zsteps
        allocate(Rg(3,npoints))
        l = 1
        do i = 1, xsteps
            do j = 1, ysteps
                do k = 1, zsteps
                    Rg(1,l) = origin(1) + (i - 1) * step(1)
                    Rg(2,l) = origin(2) + (j - 1) * step(2)
                    Rg(3,l) = origin(3) + (k - 1) * step(3)
                    l = l + 1
                end do
            end do
        end do
    end if

    call read_potential(trim(potfile))

    if (pe_mul_damp .or. pe_core_damp) then
        call setup_damping()
    end if

 10 format(4x,a)
 20 format(/6x,a)
 30 format(/8x,a)
 40 format(8x,a)
 50 format(6x,a,f8.4)
 60 format(6x,a,es8.1)

    write(luout, *)
    write(luout, 10) '.------------------------------------------------.'
    write(luout, 10)  '| Information from polarizable embedding library |'
    write(luout, 10) '`------------------------------------------------´'
    if (nsites > 0) then
        write(luout, '(/6x,a,i6)') 'Number of classical sites: ', nsites
    end if
    if (lmul(5)) then
        write(luout, 20) 'Multipole moments upto 5th order.'
    else if (lmul(4)) then
        write(luout, 20) 'Multipole moments upto 4th order.'
    else if (lmul(3)) then
        write(luout, 20) 'Multipole moments upto 3rd order.'
    else if (lmul(2)) then
        write(luout, 20) 'Multipole moments upto 2nd order.'
    else if (lmul(1)) then
        write(luout, 20) 'Multipole moments upto 1st order.'
    else if (lmul(0)) then
        write(luout, 20) 'Multipole moments upto 0th order.'
    end if
    if (pe_polar) then
        if (lpol(1,1)) then
            write(luout, 20) 'Dipole-dipole polarizabilities.'
        end if
        if (pe_gspol) then
            write(luout, 20) 'Dynamic response from environment will be&
                             & neglected during response calculation.'
        end if
        if (pe_nomb) then
            write(luout, 20) 'Many-body interactions will be neglected.'
        end if
        if (pe_iter) then
            write(luout, 20) 'Iterative solver for induced moments will&
                             & be used'
            write(luout, 60) 'with convergence threshold:', thriter
            if (pe_redthr) then
                write(luout, 20) 'Using reduced threshold in first few&
                                 & SCF iterations.'
            end if
        else
            write(luout, 20) 'Direct solver for induced moments will be used.'
        end if
        if (pe_ind_damp) then
            write(luout, 20) 'Interactions between inducible moments will&
                             & be damped using Tholes scheme.'
            write(luout, 50) 'Damping coefficient:', ind_damp
        end if
        if (pe_mul_damp) then
            write(luout, 20) 'Interactions between permanent and inducible&
                             & moments will be damped using Tholes scheme.'
            write(luout, 50) 'Damping coefficient:', mul_damp
        end if
        if (pe_core_damp) then
            write(luout, 20) 'Interactions between electrons/nuclei and&
                                   & inducible moments will be damped using&
                                   & Tholes scheme.'
            write(luout, 50) 'Damping coefficient:', core_damp
            write(luout,'(6x,a)') 'Core-region polarizabilities:'
            do i = 1, nnucs
                write(luout,'(12x,a,1x,f8.4)') charge2elem(Zm(1,i)), core_alphas(1,i)
            enddo
        endif
    end if
    if (pe_restart) then
         write(luout, 20) 'Existing files will be used to restart if possible.'
    end if
    if (pe_cube .and. .not. cube_field) then
        write(luout, 20) 'Cube files containing the potential from the&
                         & embedding potential will be written.'
    else if (pe_cube .and. cube_field) then
        write(luout, 20) 'Cube files containing the potential and electric&
                         & field from the embedding potential will be written.'
    end if

    ! handling sites near border
    ! -----------------------------------------------
    if (pe_border) then
        ! first locate all sites within given threshold of core nuclei
        allocate(idxs(nsites))
        idxs = 0; nidx = 0
        do i = 1, nnucs
            do j = 1, nsites
                lexist = .false.
                do k = 1, nidx
                    if (j == idxs(k)) then
                        lexist = .true.
                        exit
                    end if
                end do
                if (lexist) cycle
                Rsm = Rm(:,i) - Rs(:,j)
                if (nrm2(Rsm) <= Rmin) then
                    nidx = nidx + 1
                    idxs(nidx) = j
                end if
            end do
        end do

        if (border_type == 'REMOVE') then
            write(luout, *) ''
            do i = 1, nidx
                fd = '(8x,a,i6)'
                write(luout, fd) 'Removing all parameters on site:', idxs(i)
                if (lmul(0)) then
                    M0s(:,idxs(i)) = 0.0_dp
                endif
                if (lmul(1)) then
                    M1s(:,idxs(i)) = 0.0_dp
                endif
                if (lmul(2)) then
                    M2s(:,idxs(i)) = 0.0_dp
                endif
                if (lmul(3)) then
                    M3s(:,idxs(i)) = 0.0_dp
                endif
                if (lmul(4)) then
                    M4s(:,idxs(i)) = 0.0_dp
                endif
                if (lmul(5)) then
                    M5s(:,idxs(i)) = 0.0_dp
                endif
                if (lpol(1,1)) then
                    P11s(:,idxs(i)) = 0.0_dp
                end if
            end do
        else if (border_type == 'REDIST') then
            allocate(redists(nsites))
            redists = .false.
            do i = 1, nidx
                rclose = 1.0e10_dp
                do j = 1, nsites
                    lexist = .false.
                    do k = 1, nidx
                        if (j == idxs(k)) then
                            lexist = .true.
                            exit
                        end if
                    end do
                    if (lexist) cycle
                    Rji = Rs(:,idxs(i)) - Rs(:,j)
                    if (nrm2(Rji) <= rclose) then
                        rclose = nrm2(Rji)
                        idx = j
                    end if
                end do

                if (abs(redist_order) >= 1) then
                    if (lmul(0)) then
                        M0s(:,idx) = M0s(:,idx) + M0s(:,idxs(i)) / real(nredist, dp)
                    endif
                end if
                if (abs(redist_order) >= 2) then
                    if (lmul(1)) then
                        M1s(:,idx) = M1s(:,idx) + M1s(:,idxs(i)) / real(nredist, dp)
                    endif
                end if
                if (abs(redist_order) >= 3) then
                    if (lmul(2)) then
                        M2s(:,idx) = M2s(:,idx) + M2s(:,idxs(i)) / real(nredist, dp)
                    endif
                end if
                if (abs(redist_order) >= 4) then
                    if (lmul(3)) then
                        M3s(:,idx) = M3s(:,idx) + M3s(:,idxs(i)) / real(nredist, dp)
                    endif
                end if
                if (abs(redist_order) >= 5) then
                    if (lmul(4)) then
                        M4s(:,idx) = M4s(:,idx) + M4s(:,idxs(i)) / real(nredist, dp)
                    endif
                end if
                if (abs(redist_order) >= 6) then
                    if (lmul(5)) then
                        M5s(:,idx) = M5s(:,idx) + M5s(:,idxs(i)) / real(nredist, dp)
                    endif
                end if
                if (redist_order <= - 1) then
                    if (lpol(1,1)) then
                        P11s(:,idx) = P11s(:,idx) + P11s(:,idxs(i)) / real(nredist, dp)
                    end if
                end if

                if (nredist > 1) then
                    rclose = 1.0e10_dp
                    do j = 1, nsites
                        if (j == idx) cycle
                        lexist = .false.
                        do k = 1, nidx
                            if (j == idxs(k)) then
                                lexist = .true.
                                exit
                            end if
                        end do
                        if (lexist) cycle
                        Rji = Rs(:,idxs(i)) - Rs(:,j)
                        if (nrm2(Rji) <= rclose) then
                            rclose = nrm2(Rji)
                            jdx = j
                        end if
                    end do
                    if (abs(redist_order) >= 1) then
                        if (lmul(0)) then
                            M0s(:,jdx) = M0s(:,jdx) + M0s(:,idxs(i)) / real(nredist, dp)
                        endif
                    end if
                    if (abs(redist_order) >= 2) then
                        if (lmul(1)) then
                            M1s(:,jdx) = M1s(:,jdx) + M1s(:,idxs(i)) / real(nredist, dp)
                        endif
                    end if
                    if (abs(redist_order) >= 3) then
                        if (lmul(2)) then
                            M2s(:,jdx) = M2s(:,jdx) + M2s(:,idxs(i)) / real(nredist, dp)
                        endif
                    end if
                    if (abs(redist_order) >= 4) then
                        if (lmul(3)) then
                            M3s(:,jdx) = M3s(:,jdx) + M3s(:,idxs(i)) / real(nredist, dp)
                        endif
                    end if
                    if (abs(redist_order) >= 5) then
                        if (lmul(4)) then
                            M4s(:,jdx) = M4s(:,jdx) + M4s(:,idxs(i)) / real(nredist, dp)
                        endif
                    end if
                    if (abs(redist_order) >= 6) then
                        if (lmul(5)) then
                            M5s(:,jdx) = M5s(:,jdx) + M5s(:,idxs(i)) / real(nredist, dp)
                        endif
                    end if
                    if (redist_order <= - 1) then
                        if (lpol(1,1)) then
                            P11s(:,jdx) = P11s(:,jdx) + P11s(:,idxs(i)) / real(nredist, dp)
                        end if
                    end if
                end if

                if (nredist > 2) then
                    rclose = 1.0e10_dp
                    do j = 1, nsites
                        if (j == idx .or. j == jdx) cycle
                        lexist = .false.
                        do k = 1, nidx
                            if (j == idxs(k)) then
                                lexist = .true.
                                exit
                            end if
                        end do
                        if (lexist) cycle
                        Rji = Rs(:,idxs(i)) - Rs(:,j)
                        if (nrm2(Rji) <= rclose) then
                            rclose = nrm2(Rji)
                            kdx = j
                        end if
                    end do
                    if (abs(redist_order) >= 1) then
                        if (lmul(0)) then
                            M0s(:,kdx) = M0s(:,kdx) + M0s(:,idxs(i)) / real(nredist, dp)
                        endif
                    end if
                    if (abs(redist_order) >= 2) then
                        if (lmul(1)) then
                            M1s(:,kdx) = M1s(:,kdx) + M1s(:,idxs(i)) / real(nredist, dp)
                        endif
                    end if
                    if (abs(redist_order) >= 3) then
                        if (lmul(2)) then
                            M2s(:,kdx) = M2s(:,kdx) + M2s(:,idxs(i)) / real(nredist, dp)
                        endif
                    end if
                    if (abs(redist_order) >= 4) then
                        if (lmul(3)) then
                            M3s(:,kdx) = M3s(:,kdx) + M3s(:,idxs(i)) / real(nredist, dp)
                        endif
                    end if
                    if (abs(redist_order) >= 5) then
                        if (lmul(4)) then
                            M4s(:,kdx) = M4s(:,kdx) + M4s(:,idxs(i)) / real(nredist, dp)
                        endif
                    end if
                    if (abs(redist_order) >= 6) then
                        if (lmul(5)) then
                            M5s(:,kdx) = M5s(:,kdx) + M5s(:,idxs(i)) / real(nredist, dp)
                        endif
                    end if
                    if (redist_order <= - 1) then
                        if (lpol(1,1)) then
                            P11s(:,kdx) = P11s(:,kdx) + P11s(:,idxs(i)) / real(nredist, dp)
                        end if
                    end if
                end if

                if (lmul(0)) then
                    M0s(:,idxs(i)) = 0.0_dp
                endif
                if (lmul(1)) then
                    M1s(:,idxs(i)) = 0.0_dp
                endif
                if (lmul(2)) then
                    M2s(:,idxs(i)) = 0.0_dp
                endif
                if (lmul(3)) then
                    M3s(:,idxs(i)) = 0.0_dp
                endif
                if (lmul(4)) then
                    M4s(:,idxs(i)) = 0.0_dp
                endif
                if (lmul(5)) then
                    M5s(:,idxs(i)) = 0.0_dp
                endif
                if (lpol(1,1)) then
                    P11s(:,idxs(i)) = 0.0_dp
                end if

                if (redist_order == 1) then
                    write(luout, 20) 'Redistributing multipoles upto 0th order'
                else if (redist_order == 2) then
                    write(luout, 20) 'Redistributing multipoles upto 1st order'
                else if (redist_order == 3) then
                    write(luout, 20) 'Redistributing multipoles upto 2nd order'
                else if (redist_order == 4) then
                    write(luout, 20) 'Redistributing multipoles upto 3rd order'
                else if (redist_order == 5) then
                    write(luout, 20) 'Redistributing multipoles upto 4th order'
                else if (redist_order == 6) then
                    write(luout, 20) 'Redistributing multipoles upto 5th order'
                else if (redist_order == - 1) then
                    write(luout, 20) 'Redistributing multipoles upto 0th&
                                     & order and polarizabilities'
                else if (redist_order == - 2) then
                    write(luout, 20) 'Redistributing multipoles upto 1st&
                                     & order and polarizabilities'
                else if (redist_order == - 3) then
                    write(luout, 20) 'Redistributing multipoles upto 2nd&
                                     & order and polarizabilities'
                else if (redist_order == - 4) then
                    write(luout, 20) 'Redistributing multipoles upto 3rd&
                                     & order and polarizabilities'
                else if (redist_order == - 5) then
                    write(luout, 20) 'Redistributing multipoles upto 4th&
                                     & order and polarizabilities'
                else if (redist_order == - 6) then
                    write(luout, 20) 'Redistributing multipoles upto 5th&
                                     & order and polarizabilities'
                end if
                write(luout, '(4x,a,i6)') 'from site:', idxs(i)
                fd = '(8x,a,3i6)'
                if (nredist == 3) then
                    write(luout, fd) 'to neighbouring sites:', idx, jdx, kdx
                else if (nredist == 2) then
                    write(luout, fd) 'to neighbouring sites:', idx, jdx
                else if (nredist == 1) then
                    write(luout, fd) 'to neighbouring site:', idx
                end if
                redists(idx) = .true.
                if (nredist > 1) redists(jdx) = .true.
                if (nredist > 2) redists(kdx) = .true.
            end do
            if (abs(redist_order) >= 1) then
                if (lmul(0)) then
                    write(luout, 30) ' Resulting monopoles: '
                    write(luout, 40) '----------------------'
                    do i = 1, nsites
                        if (redists(i)) then
                            fd = '(8x,a,1x,i6,2x,f9.4)'
                            write(luout, fd) elems(i), i, M0s(:,i)
                        end if
                    end do
                end if
            end if
            if (abs(redist_order) >= 2) then
                if (lmul(1)) then
                    write(luout, 30) ' Resulting dipoles: '
                    write(luout, 40) '--------------------'
                    do i = 1, nsites
                        if (redists(i)) then
                            fd = '(8x,a,1x,i6,2x,3f9.4)'
                            write(luout, fd) elems(i), i, M1s(:,i)
                        end if
                    end do
                end if
            end if
            if (abs(redist_order) >= 3) then
                if (lmul(2)) then
                    write(luout, 30) ' Resulting quadrupoles: '
                    write(luout, 40) '------------------------'
                    do i = 1, nsites
                        if (redists(i)) then
                            fd = '(8x,a,1x,i6,2x,6f9.4)'
                            write(luout, fd) elems(i), i, M2s(:,i)
                        end if
                    end do
                end if
            end if
            if (abs(redist_order) >= 4) then
                if (lmul(3)) then
                    write(luout, 30) ' Resulting octopoles: '
                    write(luout, 40) '----------------------'
                    do i = 1, nsites
                        if (redists(i)) then
                            fd = '(8x,a,1x,i6,2x,10f9.4)'
                            write(luout, fd) elems(i), i, M3s(:,i)
                        end if
                    end do
                end if
            end if
            if (abs(redist_order) >= 5) then
                if (lmul(4)) then
                    write(luout, 30) ' Resulting hexadecapoles: '
                    write(luout, 40) '--------------------------'
                    do i = 1, nsites
                        if (redists(i)) then
                            fd = '(8x,a,1x,i6,2x,15f9.4)'
                            write(luout, fd) elems(i), i, M4s(:,i)
                        end if
                    end do
                end if
            end if
            if (abs(redist_order) >= 6) then
                if (lmul(5)) then
                    write(luout, 30) ' Resulting ditriacontapoles: '
                    write(luout, 40) '-----------------------------'
                    do i = 1, nsites
                        if (redists(i)) then
                            fd = '(8x,a,1x,i6,2x,21f9.4)'
                            write(luout, fd) elems(i), i, M5s(:,i)
                        end if
                    end do
                end if
            end if
            if (redist_order <= - 1) then
                if (lpol(1,1)) then
                    write(luout, 30) ' Resulting polarizabilities: '
                    write(luout, 40) '-----------------------------'
                    do i = 1, nsites
                        if (redists(i)) then
                            fd = '(8x,a,1x,i6,2x,6f9.4)'
                            write(luout, fd) elems(i), i, P11s(:,i)
                        end if
                    end do
                end if
            end if
            deallocate(redists)
        else if (border_type == 'CHGRED') then
            write(luout, *) ''
            redist_charge = 0.0_dp
            do i = 1, nidx
                fd = '(6x,a,i6)'
                if (lmul(0)) then
                    write(luout, fd) 'Redistributing charges from site:', idxs(i)
                    redist_charge = redist_charge + M0s(1,idxs(i))
                    M0s(:,idxs(i)) = 0.0_dp
                end if
                if (lmul(1)) then
                    M1s(:,idxs(i)) = 0.0_dp
                endif
                if (lmul(2)) then
                    M2s(:,idxs(i)) = 0.0_dp
                endif
                if (lmul(3)) then
                    M3s(:,idxs(i)) = 0.0_dp
                endif
                if (lmul(4)) then
                    M4s(:,idxs(i)) = 0.0_dp
                endif
                if (lmul(5)) then
                    M5s(:,idxs(i)) = 0.0_dp
                endif
                if (lpol(1,1)) then
                    P11s(:,idxs(i)) = 0.0_dp
                end if
            end do
            redist_charge = redist_charge / (nsites - nidx)
            skip = .false.
            do i = 1, nsites
                do j = 1, nidx
                    if (i == idxs(j)) then
                        skip = .true.
                        exit
                    end if
                end do
                if (skip) then
                    skip = .false.
                    cycle
                end if
                M0s(1,i) = M0s(1,i) + redist_charge
            end do
        end if
        deallocate(idxs)
        write(luout, *) ''
    end if

    if (lmul(0)) then
        redist_charge = 0.0_dp
        do i = 1, nsites
            redist_charge = redist_charge + M0s(1,i)
        end do
        fd = '(6x,a,f9.4)'
        write(luout, fd) 'Total charge of classical region: ', redist_charge
    end if

    write(luout, *) ''

    ! number of polarizabilities different from zero
    if (any(lpol)) then
        allocate(zeroalphas(nsites))
        do i = 1, nsites
            if (maxval(abs(P11s(:,i))) <= zero) then
                zeroalphas(i) = .true.
            else
                zeroalphas(i) = .false.
                npols = npols + 1
            end if
        end do
    end if

end subroutine pe_init

!------------------------------------------------------------------------------

subroutine pe_finalize()

    use pe_variables
    use pe_integral_interfaces

    initialized = .false.

    if (allocated(Rm)) deallocate(Rm)
    if (allocated(Zm)) deallocate(Zm)
    if (allocated(Rg)) deallocate(Rg)
    if (allocated(elems)) deallocate(elems)
    if (allocated(Rs)) deallocate(Rs)
    if (allocated(Zs)) deallocate(Zs)
    if (allocated(M0s)) deallocate(M0s)
    if (allocated(M1s)) deallocate(M1s)
    if (allocated(M2s)) deallocate(M2s)
    if (allocated(M3s)) deallocate(M3s)
    if (allocated(M4s)) deallocate(M4s)
    if (allocated(M5s)) deallocate(M5s)
    if (allocated(P11s)) deallocate(P11s)
    if (allocated(exclists)) deallocate(exclists)
    if (allocated(zeroalphas)) deallocate(zeroalphas)
    if (allocated(Epe)) deallocate(Epe)
    if (allocated(Ees)) deallocate(Ees)
    if (allocated(Epol)) deallocate(Epol)
    if (allocated(displs)) deallocate(displs)
    if (allocated(sitedists)) deallocate(sitedists)
    if (allocated(siteloops)) deallocate(siteloops)
    if (allocated(poldists)) deallocate(poldists)
    if (allocated(cubeloops)) deallocate(cubeloops)
    if (allocated(cubedists)) deallocate(cubedists)

    call finalize_pe_integral_interfaces

end subroutine pe_finalize

!------------------------------------------------------------------------------

subroutine pe_input_reader(word, luinp)

    use pe_variables
    use pe_constants
    use pe_utils, only: chcase

    character(len=*), intent(inout) :: word
    integer, intent(in) :: luinp

    character(len=80) :: option
    character(len=2) :: auoraa
    integer :: i, j
    real(dp), dimension(2) :: temp

    do
        read(luinp, '(a80)') option
        call chcase(option)

        ! Read potential (optionally from potfile)
        if (trim(option(2:7)) == 'POTENT') then
            read(luinp, '(a80)') option
            backspace(luinp)
            if ((option(1:1) /= '.') .and. (option(1:1) /= '*') .and.&
              & (option(1:1) /= '!') .and. (option(1:1) /= '#')) then
                read(luinp, '(a80)') potfile
            end if
        ! direct solver for induced moments
        else if (trim(option(2:7)) == 'DIRECT') then
            read(luinp, '(a80)') option
            backspace(luinp)
            if ((option(1:1) /= '.') .and. (option(1:1) /= '*') .and.&
               & (option(1:1) /= '!') .and. (option(1:1) /= '#')) then
                read(luinp, '(a80)') option
                call chcase(option)
                if (trim(option(1:7)) == 'NOCHOL') then
                    chol = .false.
                end if 
            end if
            pe_iter = .false.
        ! iterative solver for induced moments (default)
        else if (trim(option(2:7)) == 'ITERAT') then
            read(luinp, '(a80)') option
            backspace(luinp)
            if ((option(1:1) /= '.') .and. (option(1:1) /= '*') .and.&
              & (option(1:1) /= '!') .and. (option(1:1) /= '#')) then
                read(luinp, *) thriter
            end if
            pe_iter = .true.
        ! use reduced threshold in iterative induced moments solver
        else if (trim(option(2:7)) == 'REDTHR') then
            read(luinp, '(a80)') option
            backspace(luinp)
            if ((option(1:1) /= '.') .and. (option(1:1) /= '*') .and.&
               & (option(1:1) /= '!') .and. (option(1:1) /= '#')) then
                read(luinp, *) redlvl
            end if
            pe_redthr = .true.
        ! handling sites near quantum-classical border
         else if (trim(option(2:7)) == 'BORDER') then
            read(luinp, '(a80)') option
            backspace(luinp)
            if ((option(1:1) /= '.') .and. (option(1:1) /= '*') .and.&
              & (option(1:1) /= '!') .and. (option(1:1) /= '#')) then
                read(luinp, '(a)', advance='no') border_type
                backspace(luinp)
                call chcase(border_type)
                if ((trim(border_type) /= 'REMOVE') .and.&
                  & (trim(border_type) /= 'REDIST') .and.&
                  & (trim(border_type) /= 'CHGRED')) then
                    stop 'ERROR: unknown handling of border sites'
                else if (trim(border_type) == 'REMOVE') then
                    read(luinp, *) border_type, Rmin, auoraa
                else if (trim(border_type) == 'REDIST') then
                    read(luinp, *) border_type, redist_order, Rmin, auoraa, nredist
                    if ((nredist > 3) .or. (nredist < 1)) then
                        stop 'ERROR: parameters can only be distributed to&
                             & minimum one site and maximum three sites'
                    end if
                else if (trim(border_type) == 'CHGRED') then
                    read(luinp, *) border_type, Rmin, auoraa
                else
                    stop 'ERROR: unrecognized input in .BORDER option'
                end if
                call chcase(auoraa)
                if (trim(auoraa) == 'AA') Rmin = Rmin * aa2bohr
            end if
            pe_border = .true.
        ! damp electric field from induced multipoles
        else if (trim(option(2:7)) == 'DAMP I') then
            read(luinp, '(a80)') option
            backspace(luinp)
            if ((option(1:1) /= '.') .and. (option(1:1) /= '*') .and.&
              & (option(1:1) /= '!') .and. (option(1:1) /= '#')) then
                read(luinp, *) ind_damp
            end if
            pe_ind_damp = .true.
        ! damp electric field from permanent multipoles
        else if (trim(option(2:7)) == 'DAMP M') then
            read(luinp, '(a80)') option
            backspace(luinp)
            if ((option(1:1) /= '.') .and. (option(1:1) /= '*') .and.&
              & (option(1:1) /= '!') .and. (option(1:1) /= '#')) then
                read(luinp, *) mul_damp
            end if
            pe_mul_damp = .true.
        ! damp electric field from core region
        else if (trim(option(2:7)) == 'DAMP C') then
            read(luinp, '(a80)') option
            backspace(luinp)
            if ((option(1:1) /= '.') .and. (option(1:1) /= '*') .and.&
              & (option(1:1) /= '!') .and. (option(1:1) /= '#')) then
                read(luinp, *) core_damp
                ! attempt to optionally read a custom specification of the polarizabilities
                read(luinp, '(a80)') option
                backspace(luinp)
                if ((option(1:1) /= '.') .and. (option(1:1) /= '*') .and.&
                   & (option(1:1) /= '!') .and. (option(1:1) /= '#')) then
                    read(luinp, *) j
                    allocate(core_alphas(6,j))
                    core_alphas = 0.0_dp
                    do i = 1, j
                        read(luinp, *) core_alphas(1,i)
                        core_alphas(4,i) = core_alphas(1,i)
                        core_alphas(6,i) = core_alphas(1,i)
                    enddo
                end if
            end if
            pe_core_damp = .true.
        ! the old deprecated DAMP option
        else if (trim(option(2:7)) == 'DAMP') then
            read(luinp, '(a80)') option
            backspace(luinp)
            if ((option(1:1) /= '.') .and. (option(1:1) /= '*') .and.&
              & (option(1:1) /= '!') .and. (option(1:1) /= '#')) then
                read(luinp, *) ind_damp
            end if
            pe_ind_damp = .true.
            print *, 'INFO: the .DAMP option is deprecated, please use .DAMP INDUCED'
        ! neglect dynamic response from environment
        else if (trim(option(2:7)) == 'GSPOL') then
            pe_gspol = .true.
        ! neglect many-body interactions
        else if (trim(option(2:7)) == 'NOMB') then
            pe_nomb = .true.
        ! Use existing files for restart
        else if (trim(option(2:7)) == 'RESTAR') then
            pe_restart = .true.
        ! request calculation of effective dipole integrals
        else if (trim(option(2:7)) == 'EEF') then
            read(luinp, '(a80)') option
            backspace(luinp)
            if ((option(1:1) /= '.') .and. (option(1:1) /= '*') .and.&
               & (option(1:1) /= '!') .and. (option(1:1) /= '#')) then
                read(luinp, *) ncrds
                allocate(crds(3,ncrds))
                do i = 1, ncrds
                    read(luinp, *) (crds(j,i), j = 1, 3)
                end do
            end if
            pe_lf = .true.
        ! provide LJ parameters for the QM region
        else if (trim(option(2:7)) == 'LJ') then
            read(luinp, '(a80)') option
            backspace(luinp)
            if ((option(1:1) /= '.') .and. (option(1:1) /= '*') .and.&
               & (option(1:1) /= '!') .and. (option(1:1) /= '#')) then
                read(luinp, *) qmLJsites
                allocate(qmLJs(2,qmLJsites))
                do i = 1, qmLJsites
                    read(luinp, *) (temp(j), j = 1, 2)
                    qmLJs(1,i) = temp(1) * 2.0_dp
                    qmLJs(2,i) = temp(2)
                end do
                lvdw = .true.
             end if
        ! Write cube files
        else if (trim(option(2:7)) == 'CUBE') then
            read(luinp, '(a80)') option
            backspace(luinp)
            if ((option(1:1) /= '.') .and. (option(1:1) /= '*') .and.&
              & (option(1:1) /= '!') .and. (option(1:1) /= '#')) then
                do
                    read(luinp, '(a80)') option
                    call chcase(option)
                    if (trim(option(1:7)) == 'COARSE') then
                        xgrid = 3
                        ygrid = 3
                        zgrid = 3
                    else if (trim(option(1:7)) == 'MEDIUM') then
                        xgrid = 6
                        ygrid = 6
                        zgrid = 6
                    else if (trim(option(1:7)) == 'FINE') then
                        xgrid = 12
                        ygrid = 12
                        zgrid = 12
                    else if (trim(option(1:7)) == 'GRID') then
                        read(luinp, *) xsize, xgrid, ysize, ygrid, zsize, zgrid
                    else if (trim(option(1:7)) == 'FIELD') then
                        cube_field = .true.
                    else if (option(1:1) == '.' .or. option(1:1) == '*') then
                        backspace(luinp)
                        exit
                    else if (option(1:1) == '!' .or. option(1:1) == '#') then
                        cycle
                    else
                        stop 'ERROR: unknown input under .CUBE option'
                    end if
                end do
            end if
            pe_cube = .true.
        ! apply external electric field
        else if (trim(option(2:7)) == 'FIELD') then
            pe_field = .true.
            read(luinp, '(a80)') option
            backspace(luinp)
            if ((option(1:1) /= '.') .and. (option(1:1) /= '*') .and.&
               & (option(1:1) /= '!') .and. (option(1:1) /= '#')) then
                read(luinp, *) (efields(i), i = 1, 3)
            end if
        ! verbose output
        else if (trim(option(2:7)) == 'VERBOS') then
            pe_verbose = .true.
        ! debug output
        else if (trim(option(2:7)) == 'DEBUG') then
            pe_debug = .true.
            pe_verbose = .true.
        ! isotropic polarizabilities
        else if (trim(option(2:7)) == 'ISOPOL') then
            pe_isopol = .true.
        ! zero out the polarizabilities
        else if (trim(option(2:7)) == 'ZEROPO') then
            pe_zeropol = .true.
        ! zero out higher-order multipoles
        else if (trim(option(2:7)) == 'ZEROMU') then
            read(luinp, '(a80)') option
            backspace(luinp)
            if ((option(1:1) /= '.') .and. (option(1:1) /= '*') .and.&
               & (option(1:1) /= '!') .and. (option(1:1) /= '#')) then
                read(luinp, *) zeromul_order
                if (zeromul_order < 0) then
                    STOP 'ERROR: ZEROMUL order cannot be negative'
                end if
            end if
            pe_zeromul = .true.
        else if (option(1:1) == '*') then
            word = option
            exit
        else if (option(1:1) == '!' .or. option(1:1) == '#') then
            cycle
        else
            write(luout, *) 'ERROR: unknown option:', option
            stop 'ERROR: unknown option in *PEQM section'
        end if
    end do

end subroutine pe_input_reader

!------------------------------------------------------------------------------

subroutine read_potential(filename)

    use pe_variables
    use pe_constants
    use pe_utils, only: chcase, openfile, elem2charge

    character(len=*) :: filename

    integer :: i, j, s
    integer :: nlines
    integer :: lupot
    integer :: multype
    integer, dimension(2) :: poltype
    integer, dimension(:), allocatable :: itemp
    real(dp) :: trace
    real(dp), dimension(21) :: temp
    character(len=2) :: bohroraa
    character(len=80) :: word
    logical :: lexist

    inquire(file=filename, exist=lexist)
    if (lexist) then
        call openfile(filename, lupot, 'old', 'formatted')
    else
        write(luout, *) 'ERROR: potential input file not found'
        stop 'ERROR: potential input file not found'
    end if

    write(luout, '(/2x,a)') 'Reading potential input file for polarizable embedding calculation'
    print *, 'PElib: reading potential input file for polarizable embedding calculation'

    do
        read(lupot, *, end=112) word
        call chcase(word)
        if (trim(word) == '@COORDINATES') then
            read(lupot, *) nsites
            read(lupot, *) bohroraa
            allocate(elems(nsites), Zs(1,nsites), Rs(3,nsites))
            do i = 1, nsites
                read(lupot, *) elems(i), (Rs(j,i), j = 1, 3)
                Zs(1,i) = elem2charge(elems(i))
            end do
            call chcase(bohroraa)
            if (bohroraa == 'AA') then
                Rs = Rs * aa2bohr
            end if
            exit
        else
            cycle
        end if
112     stop 'ERROR: no coordinates found in potential input file'
    end do

    rewind(lupot)

    do
        read(lupot, *, end=100) word
        call chcase(word)
        if (word(1:1) /= '@') then
            cycle
        else if (trim(word) == '@COORDINATES') then
            cycle
        else if (trim(word) == '@MULTIPOLES') then
            do
                read(lupot, *, end=100) word
                call chcase(word)
                if (word(1:1) == '@') then
                    backspace(lupot)
                    exit
                else if (trim(word) /= 'ORDER') then
                    cycle
                else if (trim(word) == 'ORDER') then
                    backspace(lupot)
                    read(lupot, *) word, multype
                    if (multype == 0) then
                        if (pe_zeromul .and. (zeromul_order <= multype)) cycle
                        lmul(0) = .true.
                        allocate(M0s(1,nsites))
                        M0s = 0.0_dp
                        read(lupot, *) nlines
                        do i = 1, nlines
                            read(lupot, *) s, temp(1)
                            M0s(1,s) = temp(1)
                        end do
                    else if (multype == 1) then
                        if (pe_zeromul .and. (zeromul_order <= multype)) cycle
                        lmul(1) = .true.
                        allocate(M1s(3,nsites))
                        M1s = 0.0_dp
                        read(lupot, *) nlines
                        do i = 1, nlines
                            read(lupot, *) s, (temp(j), j = 1, 3)
                            M1s(:,s) = temp(1:3)
                        end do
                    else if (multype == 2) then
                        if (pe_zeromul .and. (zeromul_order <= multype)) cycle
                        lmul(2) = .true.
                        allocate(M2s(6,nsites))
                        M2s = 0.0_dp
                        read(lupot, *) nlines
                        do i = 1, nlines
                            read(lupot, *) s, (temp(j), j = 1, 6)
                            ! remove trace
                            trace = (temp(1) + temp(4) + temp(6)) / 3.0_dp
                            temp(1) = temp(1) - trace
                            temp(4) = temp(4) - trace
                            temp(6) = temp(6) - trace
                            M2s(:,s) = temp(1:6)
                        end do
                    else if (multype == 3) then
                        if (pe_zeromul .and. (zeromul_order <= multype)) cycle
                        lmul(3) = .true.
                        allocate(M3s(10,nsites))
                        M3s = 0.0_dp
                        read(lupot, *) nlines
                        do i = 1, nlines
                            read(lupot, *) s, (temp(j), j = 1, 10)
                            trace = (temp(1) + temp(4) + temp(6)) / 5.0_dp
                            temp(1) = temp(1) - 3.0_dp * trace
                            temp(4) = temp(4) - trace
                            temp(6) = temp(6) - trace
                            trace = (temp(2) + temp(7) + temp(9)) / 5.0_dp
                            temp(2) = temp(2) - trace
                            temp(7) = temp(7) - 3.0_dp * trace
                            temp(9) = temp(9) - trace
                            trace = (temp(3) + temp(8) + temp(10)) / 5.0_dp
                            temp(3) = temp(3) - trace
                            temp(8) = temp(8) - trace
                            temp(10) = temp(10) - 3.0_dp * trace
                            M3s(:,s) = temp(1:10)
                        end do
                    else if (multype == 4) then
                        if (pe_zeromul .and. (zeromul_order <= multype)) cycle
                        write(luout, *) 'WARNING: results will be wrong if&
                                        & non-traceless hexadecapoles&
                                        & (16-poles) are used'
                        lmul(4) = .true.
                        allocate(M4s(15,nsites))
                        M4s = 0.0_dp
                        read(lupot, *) nlines
                        do i = 1, nlines
                            read(lupot, *) s, (temp(j), j = 1, 15)
                            M4s(:,s) = temp(1:15)
                        end do
                    else if (multype == 5) then
                        if (pe_zeromul .and. (zeromul_order <= multype)) cycle
                        write(luout, *) 'WARNING: results will be wrong if&
                                        & non-traceless ditriacontapoles&
                                        & (32-poles) are used'
                        lmul(5) = .true.
                        allocate(M5s(21,nsites))
                        M5s = 0.0_dp
                        read(lupot, *) nlines
                        do i = 1, nlines
                            read(lupot, *) s, (temp(j), j = 1, 21)
                            M5s(:,s) = temp(1:21)
                        end do
                    else
                        write(luout, *) 'ERROR: unsupported or unrecognized multipole order'
                        write(luout, *) multype
                        stop 'ERROR: unsupported or unrecognized multipole order'
                    end if
                else
                    write(luout, *) 'ERROR: unknown keyword in @MULTIPOLES section'
                    write(luout, *) trim(word)
                    stop 'ERROR: unknown keyword in @MULTIPOLES section'
                end if
            end do
        else if (trim(word) == '@POLARIZABILITIES') then
            do
                read(lupot, *, end=100) word
                call chcase(word)
                if (word(1:1) == '@') then
                    backspace(lupot)
                    exit
                else if (trim(word) /= 'ORDER' .and. trim(word) /= 'EXCLISTS') then
                    cycle
                else if (trim(word) == 'ORDER') then
                    backspace(lupot)
                    read(lupot, *) word, poltype
                    if (poltype(1) == 1 .and. poltype(2) == 1) then
                        if (pe_zeropol) cycle
                        lpol(1,1) = .true.
                        pe_polar = .true.
                        if (.not. allocated(P11s)) then
                            allocate(P11s(6,nsites))
                            P11s = 0.0_dp
                        end if
                        read(lupot, *) nlines
                        do i = 1, nlines
                            read(lupot, *) s, (temp(j), j = 1, 6)
                            if (pe_isopol) then
                                trace = temp(1) + temp(4) + temp(6)
                                temp = 0.0_dp
                                temp(1) = trace / 3.0_dp
                                temp(4) = temp(1)
                                temp(6) = temp(1)
                            end if
                            P11s(:,s) = temp(1:6)
                        end do
                    else if (word(1:1) == '!' .or. word(1:1) == '#') then
                        cycle
                    else
                        write(luout, *) 'ERROR: unsupported or unrecognized polarizability order'
                        write(luout, *) poltype
                        stop 'ERROR: unsupported or unrecognized polarizability'
                    end if
                else if (trim(word) == 'EXCLISTS') then
                    if (pe_zeropol) cycle
                    read(lupot, *) nlines, lexlst
                    allocate(exclists(lexlst,nsites))
                    allocate(itemp(lexlst))
                    exclists = 0
                    do i = 1, nsites
                        exclists(1,i) = i
                    end do
                    if (lexlst > 1) then
                        do i = 1, nlines
                            itemp = 0
                            read(lupot, *) s, (itemp(j), j = 2, lexlst)
                            exclists(2:lexlst,s) = itemp(2:lexlst)
                        end do
                    end if
                    deallocate(itemp)
                else
                    write(luout, *) 'ERROR: unknown keyword in @POLARIZABILITIES section'
                    write(luout, *) trim(word)
                    stop 'ERROR: unknown keyword in @POLARIZABILITIES section'
                end if
            end do
        else if (trim(word) == '@LJ') then
            !lvdw = .true.
            allocate(LJs(2,nsites))
            LJs = 0.0_dp
            read(lupot,*) nlines
            do i = 1, nlines
                read(lupot,*) s, (temp(j), j = 1, 2)
                LJs(1,s) = temp(1) * 2.0_dp
                LJs(2,s) = temp(2)
            end do
        else
            write(luout, *) 'ERROR: unknown section in potential input file'
            write(luout, *) trim(word)
            stop 'ERROR: unknown section in potential input file'
        end if
    end do

100 continue

    close(lupot)

    ! default exclusion list (everything polarizes everything)
    if (.not. allocated(exclists)) then
        lexlst = 1
        allocate(exclists(lexlst,nsites))
        do i = 1, nsites
            exclists(1,i) = i
        end do
        if (any(lpol)) then
            write(luout, *) 'WARNING: no exclusion lists provided which means all polarizable'
            write(luout, *) '         sites are polarized by everything'
        end if
    end if

end subroutine read_potential

!------------------------------------------------------------------------------

subroutine pe_master(runtype, triang, ndim, nmats, denmats, fckmats, expvals)

    use pe_mpi
    use pe_variables
    use pe_analysis_tools
    use pe_operators

    character(*), intent(in) :: runtype
    logical, intent(in), optional :: triang
    integer, intent(in), optional :: ndim
    integer, intent(in), optional :: nmats
    real(dp), dimension(:), intent(in), optional :: denmats
    real(dp), dimension(:), intent(out), optional :: fckmats
    real(dp), dimension(:), intent(out), optional :: expvals

    ! get root id and number of MPI processes
#if defined(VAR_MPI)
    call mpi_comm_rank(comm, myid, ierr)
    call mpi_comm_size(comm, nprocs, ierr)
    master = myid
    if (master /= 0) then
        stop 'ERROR: a slave cannot be master'
    end if
#else
    myid = 0
    nprocs = 1
    master = myid
#endif

    ! in-/output consistency checks
    if ((present(denmats) .or. present(fckmats)) .and. .not. (present(ndim) .and. present(triang))) then
        stop 'ERROR: triang and ndim are required input if either denmats or fckmats are provided'
    end if

    if (present(denmats) .and. .not. present(nmats)) then
        stop 'ERROR: nmats is required input if denmats is provided'
    end if

    if (present(expvals) .and. .not. present(denmats)) then
        stop 'ERROR: expvals requires denmats'
    end if

    ! set variables according to input
    if (present(triang)) then
        trimat = triang
    else
        trimat = .true.
    end if

    if (.not. trimat) then
        stop 'ERROR: only triangular matrices in packed storage available'
    end if

    if (present(ndim)) then
        nbas = ndim
        nnbas = nbas * (nbas + 1) / 2
        n2bas = nbas * nbas
    else
        nbas = 0
        nnbas = 0
        n2bas = 0
    end if

    if (present(nmats)) then
        ndens = nmats
    else
        ndens = 1
    end if

    if (present(denmats)) then
        lden = .true.
    else
        lden = .false.
    end if

    if (present(fckmats)) then
        fckmats = 0.0_dp
        lfck = .true.
    else
        lfck = .false.
    end if

    if (present(expvals)) then
        expvals = 0.0_dp
        lexp = .true.
    else
        lexp = .false.
    end if

    ! checks according to runtype
    if (present(denmats)) then
        if (trimat) then
            if (size(denmats) < ndens * nnbas) then
                stop 'ERROR: denmats too small'
            end if
        else if (.not. trimat) then
            if (size(denmats) < ndens * n2bas) then
                stop 'ERROR: denmats too small'
            end if
        end if
    end if

    fock = .false.
    energy = .false.
    response = .false.
    london = .false.
    molgrad = .false.
    effdipole = .false.

    if (runtype == 'full_fock') then
        fock = .true.
        scf_cycle = scf_cycle + 1
        if (.not. present(fckmats)) then
            stop 'ERROR: fckmats required for runtype=full_fock'
        else if (.not. present(expvals)) then
            stop 'ERROR: expvals required for runtype=full_fock'
        end if
        if (size(fckmats) < ndens * nnbas) then
            stop 'ERROR: fckmats too small'
        else if (size(expvals) < ndens) then
            stop 'ERROR: expvals too small'
        end if
    else if (runtype == 'print_energy') then
        energy = .true.
        if (.not. present(denmats)) then
            stop 'ERROR: denmats required for runtype=print_energy'
        end if
    else if (runtype == 'dynamic_response') then
        if (pe_gspol) return
        response = .true.
        if (.not. present(denmats)) then
            stop 'ERROR: denmats required for runtype=dynamic_response'
        else if (.not. present(fckmats)) then
            stop 'ERROR: fckmats requires for runtype=dynamic_response'
        end if
        if (size(fckmats) < ndens * nnbas) then
            stop 'ERROR: fckmats too small'
        end if
    else if (runtype == 'magnetic_gradient') then
        if (ndens > 1) then
            stop 'ERROR: magnetic gradient not implemented for more than one density matrix'
        end if
        london = .true.
        if (.not. present(fckmats)) then
            stop 'ERROR: fckmats required for runtype=magnetic_gradient'
        end if
        if (size(fckmats) < 3 * nnbas) then
            stop 'ERROR: fckmats too small'
        end if
    else if (runtype == 'molecular_gradient') then
        if (ndens > 1) then
            stop 'ERROR: molecular gradient not implemented for more than one density matrix'
        end if
        molgrad = .true.
        if (present(expvals)) then
            if (size(expvals) < 3 * nnucs) then
                stop 'ERROR: expvals too small'
            end if
        end if
        if (present(fckmats)) then
            if (size(fckmats) < 3 * nnucs * nnbas) then
                stop 'ERROR: fckmats too small'
            end if
        end if
    else if (runtype == 'effdipole') then
        if (present(fckmats)) then
            if (size(fckmats) < 3 * nnbas) then
                stop 'ERROR: integrals matrix too small'
            end if
            effdipole = .true.
        end if
    else
        stop 'ERROR: unknown runtype in pe_master'
    end if

    if (nprocs == 1) then
        site_start = 1
        site_finish = nsites
        cube_start = 1
        cube_finish = npoints
    end if

#if defined(VAR_MPI)
    if (nprocs > 1) then
        call mpi_bcast(trimat, 1, lmpi, master, comm, ierr)
        call mpi_bcast(ndens, 1, impi, master, comm, ierr)
        call mpi_bcast(nbas, 1, impi, master, comm, ierr)
        call mpi_bcast(nnbas, 1, impi, master, comm, ierr)
        call mpi_bcast(n2bas, 1, impi, master, comm, ierr)
        call mpi_bcast(lden, 1, lmpi, master, comm, ierr)
        call mpi_bcast(lfck, 1, lmpi, master, comm, ierr)
        call mpi_bcast(lexp, 1, lmpi, master, comm, ierr)

        if (lden) then
            call mpi_bcast(denmats(1), nnbas * ndens, rmpi, master, comm, ierr)
        end if

        if (fock) then
            call mpi_bcast(scf_cycle, 1, impi, master, comm, ierr)
        end if

        call mpi_bcast(synced, 1, lmpi, master, comm, ierr)
        if (.not. synced) then
            call mpi_sync()
        end if
    end if
#endif

    if (fock) then
        call get_fock(denmats=denmats, fckmats=fckmats, energies=expvals)
    else if (energy) then
        call get_fock(denmats=denmats)
        if (pe_cube) then
            if (ndens > 1) then
                stop 'ERROR: CUBE not implemented for more than one density matrix (ndens>1)'
            end if
            call compute_potential_cube()
            if (cube_field) call compute_field_cube()
        end if
    else if (response) then
        call get_fock(denmats=denmats, fckmats=fckmats)
    else if (london) then
        call get_magnetic_gradient(fckmats)
    else if (molgrad) then
        if (present(expvals) .and. present(fckmats)) then
            call get_molecular_gradient(denmats=denmats, fckmats=fckmats, molgrads=expvals)
        else if (present(expvals) .and. .not. present(fckmats)) then
            call get_molecular_gradient(denmats=denmats, molgrads=expvals)
        else if (present(fckmats) .and. .not. present(expvals)) then
            call get_molecular_gradient(fckmats=fckmats)
        end if
    else if (effdipole) then
        call get_effdipole(fckmats)
    end if

end subroutine pe_master

!------------------------------------------------------------------------------

#if defined(VAR_MPI)
subroutine pe_slave(runtype)

    use pe_mpi
    use pe_variables
    use pe_analysis_tools
    use pe_operators

    character(*), intent(in) :: runtype

    real(dp), dimension(:), allocatable :: denmats
    real(dp), dimension(:), allocatable :: fckmats
    real(dp), dimension(:), allocatable :: expvals

    call mpi_comm_rank(comm, myid, ierr)
    call mpi_comm_size(comm, nprocs, ierr)
    master = 0

    fock = .false.
    energy = .false.
    response = .false.
    london = .false.
    molgrad = .false.
    effdipole = .false.

    if (runtype == 'full_fock') then
        fock = .true.
    else if (runtype == 'print_energy') then
        energy = .true.
    else if (runtype == 'dynamic_response') then
        response = .true.
    else if (runtype == 'magnetic_gradient') then
        london = .true.
    else if (runtype == 'molecular_gradient') then
        molgrad = .true.
    else if (runtype == 'effdipole') then
        effdipole = .true.
    else
        stop 'ERROR: unknown runtype in pe_slave'
    end if

    call mpi_bcast(trimat, 1, lmpi, master, comm, ierr)
    call mpi_bcast(ndens, 1, impi, master, comm, ierr)
    call mpi_bcast(nbas, 1, impi, master, comm, ierr)
    call mpi_bcast(nnbas, 1, impi, master, comm, ierr)
    call mpi_bcast(n2bas, 1, impi, master, comm, ierr)
    call mpi_bcast(lden, 1, lmpi, master, comm, ierr)
    call mpi_bcast(lfck, 1, lmpi, master, comm, ierr)
    call mpi_bcast(lexp, 1, lmpi, master, comm, ierr)

    if (lden) then
        allocate(denmats(ndens*nnbas))
        call mpi_bcast(denmats(1), nnbas * ndens, rmpi, master, comm, ierr)
    end if

    if (fock .or. response) then
        allocate(fckmats(ndens*nnbas))
    else if (london) then
        allocate(fckmats(3*n2bas))
    else if (molgrad .and. lfck) then
        allocate(fckmats(3*nnucs*nnbas))
    else if (effdipole) then
        allocate(fckmats(3*nnbas))
    end if

    if (fock) then
        allocate(expvals(ndens))
    else if (molgrad .and. lexp) then
        allocate(expvals(3*nnucs))
    end if

    if (lfck) fckmats = 0.0_dp
    if (lexp) expvals = 0.0_dp

    if (fock) then
        call mpi_bcast(scf_cycle, 1, impi, master, comm, ierr)
    end if

    call mpi_bcast(synced, 1, lmpi, master, comm, ierr)
    if (.not. synced) then
        call mpi_sync()
    end if

    if (fock) then
        call get_fock(denmats=denmats, fckmats=fckmats, energies=expvals)
        deallocate(denmats)
        deallocate(fckmats)
        deallocate(expvals)
    else if (energy) then
        call get_fock(denmats=denmats)
        if (pe_cube) then
            if (ndens > 1) then
                stop 'ERROR: CUBE not implemented for more than 1 density matrix'
            end if
            call compute_potential_cube()
            if (cube_field) call compute_field_cube()
        end if
        deallocate(denmats)
    else if (response) then
        call get_fock(denmats=denmats, fckmats=fckmats)
        deallocate(denmats)
        deallocate(fckmats)
    else if (london) then
        call get_magnetic_gradient(fckmats)
        deallocate(fckmats)
    else if (molgrad) then
        if (lexp .and. lfck) then
            call get_molecular_gradient(denmats=denmats, fckmats=fckmats, molgrads=expvals)
        else if (lexp .and. .not. lfck) then
            call get_molecular_gradient(denmats=denmats, molgrads=expvals)
        else if (lfck .and. .not. lexp) then
            call get_molecular_gradient(fckmats=fckmats)
        end if
        deallocate(denmats)
        deallocate(expvals)
    else if (effdipole) then
        call get_effdipole(fckmats)
        deallocate(fckmats)
    end if

end subroutine pe_slave

!------------------------------------------------------------------------------

subroutine mpi_sync()

    use pe_mpi
    use pe_variables
    use pe_integral_interfaces

    integer :: i, j
    integer :: quotient, remainder

    if (.not. allocated(displs)) allocate(displs(0:nprocs))
    if (.not. allocated(siteloops)) allocate(siteloops(0:nprocs))
    if (.not. allocated(sitedists)) allocate(sitedists(0:nprocs-1))
    if (myid == master) then
        quotient = nsites / nprocs
        sitedists = quotient
        if (nprocs * quotient < nsites) then
            remainder = nsites - nprocs * quotient
            do i = 1, remainder
                sitedists(i-1) = sitedists(i-1) + 1
            end do
        end if
        siteloops(0) = 0
        do i = 1, nprocs
            siteloops(i) = sum(sitedists(0:i-1))
        end do
        displs(0) = 0
        do i = 1, nprocs
            displs(i) = displs(i-1) + 3 * sitedists(i-1)
        end do
    end if

    call mpi_bcast(nsites, 1, impi, master, comm, ierr)
    call mpi_bcast(siteloops(0), nprocs + 1, impi, master, comm, ierr)
    call mpi_bcast(sitedists(0), nprocs, impi, master, comm, ierr)

    site_start = siteloops(myid) + 1
    site_finish = siteloops(myid+1)

    call mpi_bcast(nnucs, 1, impi, master, comm, ierr)

    if (myid /= master .and. .not. allocated(Zm)) allocate(Zm(1,nnucs))
    call mpi_bcast(Zm(1,1), nnucs, rmpi, master, comm, ierr)

    if (myid /= master .and. .not. allocated(Rm)) allocate(Rm(3,nnucs))
    call mpi_bcast(Rm(1,1), 3 * nnucs, rmpi, master, comm, ierr)

    if (myid /= master .and. .not. allocated(core_com)) allocate(core_com(3))
    call mpi_bcast(core_com(1), 3, rmpi, master, comm, ierr)
    if (myid /= master) call init_pe_integral_interfaces(core_com, Rm)

    if (myid /= master .and. .not. allocated(Rs)) allocate(Rs(3,nsites))
    call mpi_bcast(Rs(1,1), 3 * nsites, rmpi, master, comm, ierr)

    if (myid /= master .and. .not. allocated(Zs)) allocate(Zs(1,nsites))
    call mpi_bcast(Zs(1,1), nsites, rmpi, master, comm, ierr)

    call mpi_bcast(pe_polar, 1, lmpi, master, comm, ierr)
    call mpi_bcast(lpol(0,0), 9, lmpi, master, comm, ierr)

    if (pe_polar) then
        if (any(lpol)) then
            if (.not. allocated(poldists)) allocate(poldists(0:nprocs-1))
            if (myid == master) then
                poldists = 0
                do i = 1, nprocs
                    do j = siteloops(i-1) + 1, siteloops(i)
                        if (zeroalphas(j)) then
                            continue
                        else
                            poldists(i-1) = poldists(i-1) + 1
                        end if
                    end do
                end do
            end if
            poldists = 3 * poldists
            call mpi_bcast(poldists(0), nprocs, impi, master, comm, ierr)
            call mpi_bcast(npols, 1, impi, master, comm, ierr)
            call mpi_bcast(lexlst, 1, impi, master, comm, ierr)
            if (myid /= master .and. .not. allocated(exclists)) then
                allocate(exclists(lexlst,nsites))
            end if
            call mpi_bcast(exclists(1,1), lexlst * nsites, impi, master, comm, ierr)
            if (myid /= master .and. .not. allocated(zeroalphas)) then
                allocate(zeroalphas(nsites))
            end if
            call mpi_bcast(zeroalphas(1), nsites, lmpi, master, comm, ierr)
        end if
        call mpi_bcast(pe_iter, 1, lmpi, master, comm, ierr)
        if (pe_iter .and. any(lpol)) then
            if (myid /= master .and. .not. allocated(P11s)) allocate(P11s(6,nsites))
            call mpi_bcast(P11s(1,1), 6 * nsites, rmpi, master, comm, ierr)
        end if
        call mpi_bcast(pe_nomb, 1, lmpi, master, comm, ierr)
        call mpi_bcast(pe_ind_damp, 1, lmpi, master, comm, ierr)
        call mpi_bcast(pe_mul_damp, 1, lmpi, master, comm, ierr)
        call mpi_bcast(pe_core_damp, 1, lmpi, master, comm, ierr)
        call mpi_bcast(ind_damp, 1, rmpi, master, comm, ierr)
        call mpi_bcast(mul_damp, 1, rmpi, master, comm, ierr)
        call mpi_bcast(core_damp, 1, rmpi, master, comm, ierr)
        if (pe_core_damp) then
            if (myid /= master .and. .not. allocated(core_alphas)) allocate(core_alphas(6,nnucs))
            call mpi_bcast(core_alphas, 6 * nnucs, rmpi, master, comm, ierr)
        endif
    end if

    call mpi_bcast(pe_restart, 1, lmpi, master, comm, ierr)

    call mpi_bcast(pe_lf, 1, lmpi, master, comm, ierr)
    if (pe_lf) then
        call mpi_bcast(ncrds, 1, impi, master, comm, ierr)
        if (myid /= master .and. (ncrds > 0) .and. .not. allocated(crds)) allocate(crds(3,ncrds))
        call mpi_bcast(crds(1,1), 3 * ncrds, rmpi, master, comm, ierr)
    end if

    call mpi_bcast(lmul(0), 6, lmpi, master, comm, ierr)

    if (lmul(0)) then
        if (myid /= master .and. .not. allocated(M0s)) allocate(M0s(1,nsites))
        call mpi_bcast(M0s(1,1), nsites, rmpi, master, comm, ierr)
    end if
    if (lmul(1)) then
        if (myid /= master .and. .not. allocated(M1s)) allocate(M1s(3,nsites))
        call mpi_bcast(M1s(1,1), 3 * nsites, rmpi, master, comm, ierr)
    end if
    if (lmul(2)) then
        if (myid /= master .and. .not. allocated(M2s)) allocate(M2s(6,nsites))
        call mpi_bcast(M2s(1,1), 6 * nsites, rmpi, master, comm, ierr)
    end if
    if (lmul(3)) then
        if (myid /= master .and. .not. allocated(M3s)) allocate(M3s(10,nsites))
        call mpi_bcast(M3s(1,1), 10 * nsites, rmpi, master, comm, ierr)
    end if
    if (lmul(4)) then
        if (myid /= master .and. .not. allocated(M4s)) allocate(M4s(15,nsites))
        call mpi_bcast(M4s(1,1), 15 * nsites, rmpi, master, comm, ierr)
    end if
    if (lmul(5)) then
        if (myid /= master .and. .not. allocated(M5s)) allocate(M5s(21,nsites))
        call mpi_bcast(M5s(1,1), 21 * nsites, rmpi, master, comm, ierr)
    end if

    call mpi_bcast(lvdw, 1, lmpi, master, comm, ierr)
    if (lvdw) then
        if (myid /= master .and. .not. allocated(LJs)) allocate(LJs(2,nsites))
        call mpi_bcast(LJs(1,1), 2 * nsites, rmpi, master, comm, ierr)
        call mpi_bcast(qmLJsites, 1, impi, master, comm, ierr)
        if (myid /= master .and. .not. allocated(qmLJs)) allocate(qmLJs(2,qmLJsites))
        call mpi_bcast(qmLJs(1,1), 2 * qmLJsites, rmpi, master, comm, ierr)
    end if

    call mpi_bcast(pe_cube, 1, lmpi, master, comm, ierr)

    if (pe_cube) then
        call mpi_bcast(cube_field, 1, lmpi, master, comm, ierr)
        if (.not. allocated(cubeloops)) allocate(cubeloops(0:nprocs))
        if (.not. allocated(cubedists)) allocate(cubedists(0:nprocs-1))
        if (myid == master) then
            quotient = npoints / nprocs
            cubedists = quotient
            if (nprocs * quotient < npoints) then
                remainder = npoints - nprocs * quotient
                do i = 1, remainder
                    cubedists(i-1) = cubedists(i-1) + 1
                end do
            end if
            cubeloops(0) = 0
            do i = 1, nprocs
                cubeloops(i) = sum(cubedists(0:i-1))
            end do
        end if
        call mpi_bcast(npoints, 1, impi, master, comm, ierr)
        call mpi_bcast(cubedists(0), nprocs, impi, master, comm, ierr)
        call mpi_bcast(cubeloops(0), nprocs + 1, impi, master, comm, ierr)
        cube_start = cubeloops(myid) + 1
        cube_finish = cubeloops(myid+1)
        if (myid /= master .and. .not. allocated(Rg)) allocate(Rg(3,npoints))
        call mpi_bcast(Rg(1,1), 3 * npoints, rmpi, master, comm, ierr)
    end if

    synced = .true.

end subroutine mpi_sync
#endif

!------------------------------------------------------------------------------

subroutine setup_damping()

    use pe_variables

    integer :: i

    ! check if there are polarizabilities on all sites needed for damping field
    ! from permanent multipoles
    if (pe_mul_damp) then
        do i = 1, nsites
            if (maxval(abs(P11s(:,i))) <= zero) then
                stop 'ERROR: damping field from permanent multipoles requires polarizabilities on all sites'
            end if
        end do
    end if

    ! if core polarizabilities were not read in, then attempt to use values from
    ! table 7 in Piet Th. van Duijnen and Marcel Swart, JPC A 102, 2399, 1998
    ! DOI: 10.1021/jp980221f
    if (pe_core_damp .and. .not. allocated(core_alphas)) then
        allocate(core_alphas(6,nnucs))
        core_alphas = 0.0_dp
        do i = 1, nnucs
            if (nint(Zm(1,i)) == 1) then
                core_alphas(1,i) = 2.7927_dp
            end if
            if (nint(Zm(1,i)) == 6) then
                core_alphas(1,i) = 8.6959_dp
            end if
            if (nint(Zm(1,i)) == 7) then
                core_alphas(1,i) = 6.5565_dp
            end if
            if (nint(Zm(1,i)) == 8) then
                core_alphas(1,i) = 5.7494_dp
            end if
            if (nint(Zm(1,i)) == 9) then
                core_alphas(1,i) = 3.0013_dp
            end if
            if (nint(Zm(1,i)) == 16) then
                core_alphas(1,i) = 16.6984_dp
            end if
            if (nint(Zm(1,i)) == 17) then
                core_alphas(1,i) = 16.1979_dp
            end if
            if (nint(Zm(1,i)) == 35) then
                core_alphas(1,i) = 23.5714_dp
            end if
            if (nint(Zm(1,i)) == 53) then
                core_alphas(1,i) = 36.9880_dp
            end if
            if (core_alphas(1,i) <= zero) then
                write(luout, '(a,i5)') 'ERROR: isotropic polarizability not available for core atom ', i
                stop 'ERROR: damping electric field from core failed'
            end if
            core_alphas(4,i) = core_alphas(1,i)
            core_alphas(6,i) = core_alphas(1,i)
        enddo
    else if (pe_core_damp .and. allocated(core_alphas)) then
        if (size(core_alphas, 2) /= nnucs) then
            stop 'ERROR: number of core polarizabilities given in input does not match number of core atoms'
        end if
    end if

end subroutine setup_damping

!------------------------------------------------------------------------------

end module polarizable_embedding
