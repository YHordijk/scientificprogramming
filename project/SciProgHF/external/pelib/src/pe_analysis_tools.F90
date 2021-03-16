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
module pe_analysis_tools

    use pe_precision

    implicit none

    private

    public :: compute_potential_cube, compute_field_cube

contains

!------------------------------------------------------------------------------

subroutine compute_potential_cube()

    use pe_mpi
    use pe_variables
    use pe_utils, only: openfile
    use pe_multipole_tools, only: multipole_derivative
    use pe_induced_moments, only: induced_moments
    use pe_potential_derivatives

    integer :: i, j, k, l
    integer :: lu
    logical :: lexist
    character(len=22) :: filename
    real(dp), dimension(1) :: Vp
    real(dp), dimension(3) :: Rsg
    real(dp), dimension(:), allocatable :: Vpe
    real(dp), dimension(:,:), allocatable :: M1inds
    real(dp), dimension(:), allocatable :: Fnucs, Fmuls
    real(dp), dimension(:,:), allocatable :: Fels

    if (myid == master) then
        allocate(Vpe(npoints))
    else if (myid /= master) then
        allocate(Vpe(cubedists(myid)))
    end if

    Vpe = 0.0_dp

    if (any(lmul)) then
        k = 1
        do i = cube_start, cube_finish
            Vp = 0.0_dp
            do j = 1, nsites
                Rsg = Rg(:,i) - Rs(:,j)
                if (lmul(0)) then
                    call multipole_derivative(Vp, Rsg, M0s(:,j))
                end if
                if (lmul(1)) then
                    call multipole_derivative(Vp, Rsg, M1s(:,j))
                end if
                if (lmul(2)) then
                    call multipole_derivative(Vp, Rsg, M2s(:,j))
                end if
                if (lmul(3)) then
                    call multipole_derivative(Vp, Rsg, M3s(:,j))
                end if
                if (lmul(4)) then
                    call multipole_derivative(Vp, Rsg, M4s(:,j))
                end if
                if (lmul(5)) then
                    call multipole_derivative(Vp, Rsg, M5s(:,j))
                end if
            end do
            Vpe(k) = Vpe(k) + Vp(1)
            k = k + 1
        end do
    end if

    if (lpol(1,1)) then
        allocate(M1inds(3*npols,1))
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
                read(lu) M1inds
                close(lu)
            end if
        else
            write(luout, *) 'WARNING: cannot create cube file because pe_induced_moments.bin not found'
            return
        end if
#if defined(VAR_MPI)
        if (nprocs > 1) then
            call mpi_bcast(M1inds(1,1), 3 * npols, rmpi, master, comm, ierr)
        end if
#endif
        k = 1
        do i = cube_start, cube_finish
            l = 1
            Vp = 0.0_dp
            do j = 1, nsites
                if (zeroalphas(j)) cycle
                Rsg = Rg(:,i) - Rs(:,j)
                call multipole_derivative(Vp, Rsg, M1inds(l:l+2,1))
                l = l + 3
            end do
            Vpe(k) = Vpe(k) + Vp(1)
            k = k + 1
        end do
        deallocate(M1inds)
    end if

#if defined(VAR_MPI)
    if (myid == master .and. nprocs > 1) then
        displs(0) = 0
        do i = 1, nprocs
            displs(i) = displs(i-1) + cubedists(i-1)
        end do
        call mpi_gatherv(mpi_in_place, dummy_int, rmpi, Vpe, cubedists, displs, rmpi, master, comm, ierr)
    else if (myid /= master) then
        call mpi_gatherv(Vpe(1), cubedists(myid), rmpi, dummy_real, cubedists, displs, rmpi, master, comm, ierr)
    end if
#endif
    if (myid == master) then
        call openfile('embedding_potential.cube', lu, 'unknown', 'formatted')
        write(lu, '(a)') 'Embedding potential'
        write(lu, '(a)') 'Generated by the PE library'
        write(lu, '(i5,3f12.6)') nnucs, origin
        write(lu, '(i5,3f12.6)') xsteps, step(1), 0.0_dp, 0.0_dp
        write(lu, '(i5,3f12.6)') ysteps, 0.0_dp, step(2), 0.0_dp
        write(lu, '(i5,3f12.6)') zsteps, 0.0_dp, 0.0_dp, step(3)
        do j = 1, nnucs
            write(lu, '(i5,4f12.6)') int(Zm(1,j)), Zm(1,j), Rm(:,j)
        end do
        do i = 1, xsteps * ysteps
            j = (i - 1) * zsteps + 1
            k = j - 1 + zsteps
            write(lu, '(6e13.5)') Vpe(j:k)
        end do
        close(lu)
    end if

    deallocate(Vpe)

end subroutine compute_potential_cube

!------------------------------------------------------------------------------

subroutine compute_field_cube()

    use pe_mpi
    use pe_variables
    use pe_utils, only: openfile
    use pe_multipole_tools, only: multipole_derivative

    character(len=1) :: tcl
    character(len=99) :: cl
    integer :: i, j, k, l
    integer :: lu
    logical :: lexist
    real(dp), dimension(3) :: Rsg, Fp
    real(dp), dimension(:,:), allocatable :: Fpe, M1inds


    if (myid == master) then
        allocate(Fpe(npoints,3))
    else if (myid /= master) then
        allocate(Fpe(cubedists(myid),3))
    end if

    Fpe = 0.0_dp

    if (any(lmul)) then
        k = 1
        do i = cube_start, cube_finish
            do j = 1, nsites
                Rsg = Rg(:,i) - Rs(:,j)
                if (lmul(0)) then
                    Fp = 0.0_dp
                    call multipole_derivative(Fp, Rsg, M0s(:,j))
                    Fpe(k,:) = Fpe(k,:) + Fp
                end if
                if (lmul(1)) then
                    Fp = 0.0_dp
                    call multipole_derivative(Fp, Rsg, M1s(:,j))
                    Fpe(k,:) = Fpe(k,:) + Fp
                end if
                if (lmul(2)) then
                    Fp = 0.0_dp
                    call multipole_derivative(Fp, Rsg, M2s(:,j))
                    Fpe(k,:) = Fpe(k,:) + Fp
                end if
                if (lmul(3)) then
                    Fp = 0.0_dp
                    call multipole_derivative(Fp, Rsg, M3s(:,j))
                    Fpe(k,:) = Fpe(k,:) + Fp
                end if
                if (lmul(4)) then
                    Fp = 0.0_dp
                    call multipole_derivative(Fp, Rsg, M4s(:,j))
                    Fpe(k,:) = Fpe(k,:) + Fp
                end if
                if (lmul(5)) then
                    Fp = 0.0_dp
                    call multipole_derivative(Fp, Rsg, M5s(:,j))
                    Fpe(k,:) = Fpe(k,:) + Fp
                end if
            end do
            k = k + 1
        end do
    end if

    if (lpol(1,1)) then
        allocate(M1inds(3*npols,1))
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
                read(lu) M1inds
                close(lu)
            end if
        else
            write(luout, *) 'WARNING: cannot create cube file because pe_induced_moments.bin not found'
            return
        end if
#if defined(VAR_MPI)
        if (nprocs > 1) then
            call mpi_bcast(M1inds(1,1), 3 * npols, rmpi, master, comm, ierr)
        end if
#endif
        k = 1
        do i = cube_start, cube_finish
            l = 1
            do j = 1, nsites
                if (zeroalphas(j)) cycle
                Rsg = Rg(:,i) - Rs(:,j)
                Fp = 0.0_dp
                call multipole_derivative(Fp, Rsg, M1inds(l:l+2,1))
                Fpe(k,:) = Fpe(k,:) + Fp
                l = l + 3
            end do
            k = k + 1
        end do
        deallocate(M1inds)
    end if

#if defined(VAR_MPI)
    if (myid == master .and. nprocs > 1) then
        displs(0) = 0
        do i = 1, nprocs
            displs(i) = displs(i-1) + cubedists(i-1)
        end do
        do i = 1, 3
            call mpi_gatherv(mpi_in_place, dummy_int, rmpi, Fpe(1,i), cubedists, displs, rmpi, master, comm, ierr)
        end do
    else if (myid /= master) then
        do i = 1, 3
            call mpi_gatherv(Fpe(1,i), cubedists(myid), rmpi, dummy_real, cubedists, displs, rmpi, master, comm, ierr)
        end do
    end if
#endif
    if (myid == master) then
        do l = 1, 3
            write(cl, *) l
            tcl = trim(adjustl(cl))
            call openfile('embedding_field_'//tcl//'.cube', lu, 'unknown', 'formatted')
            write(lu, '(a)') 'Embedding potential electric field component '//tcl
            write(lu, '(a)') 'Generated by the PE library'
            write(lu, '(i5,3f12.6)') nnucs, origin
            write(lu, '(i5,3f12.6)') xsteps, step(1), 0.0_dp, 0.0_dp
            write(lu, '(i5,3f12.6)') ysteps, 0.0_dp, step(2), 0.0_dp
            write(lu, '(i5,3f12.6)') zsteps, 0.0_dp, 0.0_dp, step(3)
            do j = 1, nnucs
                write(lu, '(i5,4f12.6)') nint(Zm(1,j)), Zm(1,j), Rm(:,j)
            end do
            do i = 1, xsteps * ysteps
                j = (i - 1) * zsteps + 1
                k = j - 1 + zsteps
                write(lu, '(6e13.5)') Fpe(j:k,l)
            end do
            close(lu)
        end do
    end if

    deallocate(Fpe)

end subroutine compute_field_cube

!------------------------------------------------------------------------------

end module pe_analysis_tools
