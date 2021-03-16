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
module pe_multipole_tools

    use pe_precision

    implicit none

    private

    public :: multipole_derivative, multipole_derivative_damped
    public :: T, Td, Tk_tensor, Tk_damped_tensor, prefactors, damping_coefficient

    ! C. E. Dykstra, J. Comp. Chem., 9 (1988), 476
    ! C^(n)_ij coefficients for calculating T(k) tensor elements
    integer, dimension(:,:,:), allocatable, save :: Cnij

contains

!------------------------------------------------------------------------------

subroutine multipole_derivative(Fi, Rji, Mkj)

    real(dp), dimension(:), intent(inout) :: Fi
    real(dp), dimension(3), intent(in) :: Rji
    real(dp), dimension(:), intent(in) :: Mkj

    integer :: i, j, k, l, m
    integer :: a, b, c, x, y, z
    real(dp) :: taylor, symfac
    real(dp), dimension(:), allocatable :: Tk

    ! the order of the incoming multipole moment
    k = int(0.5_dp * (sqrt(1.0_dp + 8.0_dp * real(size(Mkj), dp)) - 1.0_dp)) - 1

    ! the order of the derivative
    l = int(0.5_dp * (sqrt(1.0_dp + 8.0_dp * real(size(Fi), dp)) - 1.0_dp)) - 1

    if (mod(k + l, 2) == 0) then
        taylor = 1.0_dp / real(factorial(k), dp)
    else if (mod(k + l, 2) /= 0) then
        taylor = - 1.0_dp / real(factorial(k), dp)
    end if

    allocate(Tk((k + l + 1) * (k + l + 2) / 2))
    call Tk_tensor(Tk, Rji)

    do x = k + l, 0, - 1
        do y = k + l, 0, - 1
            do z = k + l, 0, - 1
                if (x + y + z /= k + l) cycle
                i = xyz2idx(x, y, z)
                do a = x, 0, -1
                    do b = y, 0, -1
                        do c = z, 0, -1
                            if (a + b + c /= k) cycle
                            j = xyz2idx(a, b, c)
                            m = xyz2idx(x-a, y-b, z-c)
                            symfac = real(trinom(a, b, c), dp)
                            Fi(m) = Fi(m) + taylor * symfac * Tk(i) * Mkj(j)
                        end do
                    end do
                end do
            end do
        end do
    end do

end subroutine multipole_derivative

!------------------------------------------------------------------------------

subroutine multipole_derivative_damped(Fi, Rji, Mkj, alpha_i, alpha_j, factor)

    real(dp), dimension(3), intent(inout) :: Fi
    real(dp), dimension(3), intent(in) :: Rji
    real(dp), dimension(:), intent(in) :: Mkj
    real(dp), intent(in) :: alpha_i, alpha_j, factor

    integer :: i, j, k, l, m
    integer :: a, b, c, x, y, z
    real(dp) :: taylor, symfac
    real(dp), dimension(:), allocatable :: Tk

    ! the order of the incoming multipole moment
    k = int(0.5_dp * (sqrt(1.0_dp + 8.0_dp * real(size(Mkj), dp)) - 1.0_dp)) - 1

    ! the order of the derivative
    l = int(0.5_dp * (sqrt(1.0_dp + 8.0_dp * real(size(Fi), dp)) - 1.0_dp)) - 1

    if (mod(k + l, 2) == 0) then
        taylor = 1.0_dp / real(factorial(k), dp)
    else if (mod(k + l, 2) /= 0) then
        taylor = - 1.0_dp / real(factorial(k), dp)
    end if

    allocate(Tk((k + l + 1) * (k + l + 2) / 2))
    call Tk_damped_tensor(Tk, Rji, alpha_i, alpha_j, factor)

    do x = k + l, 0, - 1
        do y = k + l, 0, - 1
            do z = k + l, 0, - 1
                if (x + y + z /= k + l) cycle
                i = xyz2idx(x, y, z)
                do a = x, 0, -1
                    do b = y, 0, -1
                        do c = z, 0, -1
                            if (a + b + c /= k) cycle
                            j = xyz2idx(a, b, c)
                            m = xyz2idx(x-a, y-b, z-c)
                            symfac = real(trinom(a, b, c), dp)
                            Fi(m) = Fi(m) + taylor * symfac * Tk(i) * Mkj(j)
                        end do
                    end do
                end do
            end do
        end do
    end do

end subroutine multipole_derivative_damped

!------------------------------------------------------------------------------

subroutine Tk_tensor(Tk, Rij)

    real(dp), dimension(:), intent(out) :: Tk
    real(dp), dimension(3), intent(in) :: Rij

    integer :: k, i
    integer :: x, y, z

    k = int(0.5_dp * (sqrt(1.0_dp + 8.0_dp * real(size(Tk), dp)) - 1.0_dp)) - 1

    i = 1
    do x = k, 0, - 1
        do y = k, 0, - 1
            do z = k, 0, - 1
                if (x + y + z /= k) cycle
                Tk(i) = T(Rij, x, y, z)
                i = i + 1
            end do
        end do
    end do

end subroutine Tk_tensor

!------------------------------------------------------------------------------

subroutine Tk_damped_tensor(Tk, Rij, alpha_i, alpha_j, factor)

    real(dp), dimension(:), intent(out) :: Tk
    real(dp), dimension(3), intent(in) :: Rij
    real(dp), intent(in) :: alpha_i, alpha_j, factor

    integer :: k, i
    integer :: x, y, z

    k = int(0.5_dp * (sqrt(1.0_dp + 8.0_dp * real(size(Tk), dp)) - 1.0_dp)) - 1

    i = 1
    do x = k, 0, - 1
        do y = k, 0, - 1
            do z = k, 0, - 1
                if (x + y + z /= k) cycle
                Tk(i) = Td(Rij, x, y, z, alpha_i, alpha_j, factor)
                i = i + 1
            end do
        end do
    end do

end subroutine Tk_damped_tensor

!------------------------------------------------------------------------------

function T(Rij, x, y, z)

    ! C. E. Dykstra, J. Comp. Chem., 9 (1988), 476

    use pe_blas_interfaces, only: nrm2

    integer, intent(in) :: x, y, z
    real(dp), dimension(3), intent(in) :: Rij

    integer :: l, m, n
    real(dp) :: T
    real(dp) :: R, Cx, Cy, Cz

    if (.not. allocated(Cnij)) call Tk_coefficients

    R = nrm2(Rij)

    T = 0.0_dp

    do l = 0, x
        Cx = Cnij(1,x,l) * (Rij(1) / R)**l
        do m = 0, y
            Cy = Cx * Cnij(l+x+1,y,m) * (Rij(2) / R)**m
            do n = 0, z
                Cz = Cy * Cnij(l+x+m+y+1,z,n) * (Rij(3) / R)**n
                T = T + Cz
            end do
        end do
    end do

    T = T / R**(x + y + z + 1)

end function T

!------------------------------------------------------------------------------

function Td(Rij, x, y, z, alpha_i, alpha_j, factor)

    ! C. E. Dykstra, J. Comp. Chem., 9 (1988), 476

    use pe_blas_interfaces, only: nrm2

    integer, intent(in) :: x, y, z
    real(dp), intent(in) :: alpha_i, alpha_j, factor
    real(dp), dimension(3), intent(in) :: Rij

    integer :: l, m, n, k
    real(dp) :: Td
    real(dp) :: R, Cx, Cy, Cz
    real(dp) :: v, fV, fE, fT, fD

    k = x + y + z

    R = nrm2(Rij)

    if (.not. allocated(Cnij)) call Tk_coefficients

    ! Thole damping
    ! JPC A 102 (1998) 2399 and Mol. Sim. 32 (2006) 471
    ! v = factor * u, where factor = 2.1304 (default) and u = R / (alpha_i * alpha_j)**(1/6)
    call damping_coefficient(Rij, alpha_i, alpha_j, factor, v)
    if (k >= 0) then
        fV = 1.0_dp - (0.5_dp * v + 1.0_dp) * exp(-v)
    end if
    if (k >= 1) then
        fE = fV  - (0.5_dp * v**2 + 0.5_dp * v) * exp(-v)
    end if
    if (k >= 2) then
        fT = fE - (1.0_dp / 6.0_dp) * v**3 * exp(-v)
    end if
    if (k >= 3) then
        fD = fT - (1.0_dp / 30.0_dp) * v**4 * exp(-v)
    end if
    if (k > 3) then
        stop 'ERROR: damping only works upto third-order interaction tensors'
    end if

    Td = 0.0_dp

    do l = 0, x
        Cx = Cnij(1,x,l) * (Rij(1) / R)**l
        do m = 0, y
            Cy = Cx * Cnij(l+x+1,y,m) * (Rij(2) / R)**m
            do n = 0, z
                Cz = Cy * Cnij(l+x+m+y+1,z,n) * (Rij(3) / R)**n
                if (l + m + n == 0) then
                    if (k == 0) then
                        Cz = Cz * fV
                    else if (k == 2) then
                        Cz = Cz * fE
                    end if
                else if (l + m + n == 1) then
                    if (k == 1) then
                        Cz = Cz * fE
                    else if (k == 3) then
                        Cz = Cz * fD
                    end if
                else if (l + m + n == 2) then
                    if (k == 2) then
                        Cz = Cz * fT
                    end if
                else if (l + m + n == 3) then
                    if (k == 3) then
                        Cz = Cz * fD
                    end if
                end if
                Td = Td + Cz
            end do
        end do
    end do

    Td = Td / R**(x + y + z + 1)

end function Td

!------------------------------------------------------------------------------

subroutine Tk_coefficients

    ! C. E. Dykstra, J. Comp. Chem., 9 (1988), 476

    integer :: i, j, k, n

! TODO
!    i = max(mulorder, polorder)
!    allocate(Cnij(2*i+3,0:i+1,0:i+1))
    allocate(Cnij(2*5+3,0:5+1,0:5+1))

    Cnij = 0
    Cnij(:,0,0) = 1
    do n = 1, 2*5+3
        if (mod(n,2) == 0) cycle
        do i = 1, 5+1
            if (mod(i,2) /= 0) then
                k = i - 1
            else if (mod(i,2) == 0) then
                k = i
            end if
            do j = 0, i
                if (mod(i+j,2) /= 0) cycle
                if (j == 0) then
                    Cnij(n,i,j) = Cnij(n,i-1,j+1)
                else if (j /= i) then
                    Cnij(n,i,j) = (j + 1) * Cnij(n,i-1,j+1)
                    Cnij(n,i,j) = Cnij(n,i,j) - (n + k) * Cnij(n,i-1,j-1)
                    k = k + 2
                else if (j == i) then
                    Cnij(n,i,j) = - (n + k) * Cnij(n,i-1,j-1)
                end if
            end do
        end do
    end do

end subroutine Tk_coefficients

!------------------------------------------------------------------------------

function xyz2idx(x, y, z) result(idx)

    integer, intent(in) :: x, y, z

    integer :: k, a, b, c
    integer :: idx

    k = x + y + z
    idx = 1
    do a = k, 0, -1
        do b = k, 0, -1
            do c = k, 0, -1
                if (a + b + c /= k) cycle
                if (a /= x .or. b /= y .or. c /= z) then
                    idx = idx + 1
                else
                    return
                end if
            end do
        end do
    end do

end function xyz2idx

!------------------------------------------------------------------------------

subroutine prefactors(factors)

    real(dp), dimension(:), intent(out) :: factors

    integer :: k
    real(dp) :: taylor

    k = int(0.5_dp * (sqrt(1.0_dp + 8.0_dp * real(size(factors), dp)) - 1.0_dp)) - 1

    if (mod(k,2) == 0) then
        taylor = 1.0_dp / real(factorial(k), dp)
    else if (mod(k,2) /= 0) then
        taylor = - 1.0_dp / real(factorial(k), dp)
    end if

    call symmetry_factors(factors)

    factors = taylor * factors

end subroutine prefactors

!------------------------------------------------------------------------------

subroutine symmetry_factors(factors)

    real(dp), dimension(:), intent(out) :: factors

    integer :: idx, x, y, z, k

    k = int(0.5_dp * (sqrt(1.0_dp + 8.0_dp * real(size(factors), dp)) - 1.0_dp)) - 1

    idx = 1
    do x = k, 0, - 1
        do y = k, 0, - 1
            do z = k, 0, - 1
                if (x + y + z /= k) cycle
                factors(idx) = real(trinom(x, y, z), dp)
                idx = idx + 1
            end do
        end do
     end do

end subroutine symmetry_factors

!------------------------------------------------------------------------------

function trinom(i, j, k)

    ! trinomial coefficient

    integer, intent(in) :: i, j, k

    integer :: trinom

    trinom = factorial(i+j+k) / (factorial(i) * factorial(j) * factorial(k))

end function trinom

!------------------------------------------------------------------------------

subroutine damping_coefficient(Rij, alpha_i, alpha_j, factor, coeff)

    use pe_variables
    use pe_blas_interfaces, only: nrm2

    real(dp), dimension(3), intent(in) :: Rij
    real(dp), intent(in) :: alpha_i, alpha_j
    real(dp), intent(in) :: factor
    real(dp), intent(out) :: coeff

    real(dp), parameter :: i6 = 1.0_dp / 6.0_dp

    ! Thole damping
    ! JPC A 102 (1998) 2399 and Mol. Sim. 32 (2006) 471
    ! coeff = a * u , where a = 2.1304 (default) and u = R / (alpha_i * alpha_j)**(1/6)

    if ((alpha_i < zero) .or. (alpha_j < zero)) then
        stop 'ERROR: damping electric fields requires non-zero polarizabilities'
    end if
    coeff = factor * nrm2(Rij) * (alpha_i * alpha_j)**(- i6)

end subroutine damping_coefficient

!------------------------------------------------------------------------------

recursive function factorial(n) result(nfact)

    ! Clive Page, http://www.star.le.ac.uk/~cgp/f90course/f90.html

    integer, intent(in) :: n

    integer :: nfact

    if (n > 0) then
        nfact = n * factorial(n-1)
    else
        nfact = 1
    end if

end function factorial

!------------------------------------------------------------------------------

end module pe_multipole_tools
