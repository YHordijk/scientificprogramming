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

module pe_blas_interfaces

    use pe_precision

    implicit none

contains

!------------------------------------------------------------------------------

function nrm2(x)

    real(dp), external :: dnrm2

    real(dp) :: nrm2
    real(dp), dimension(:), intent(in) :: x

    integer :: n, incx

    incx = 1

    n = size(x)

    nrm2 = dnrm2(n, x(1), incx)

end function nrm2

!------------------------------------------------------------------------------

function dot(x,y)

    real(dp), external :: ddot

    real(dp) :: dot
    real(dp), dimension(:), intent(in) :: x, y

    integer :: n, incx, incy

    incx = 1
    incy = 1

    n = size(x)

    dot = ddot(n, x(1), incx, y(1), incy)

end function dot

!------------------------------------------------------------------------------

subroutine axpy(x, y, a)

    external :: daxpy

    real(dp), intent(in), optional :: a
    real(dp), dimension(:), intent(in) :: x
    real(dp), dimension(:), intent(inout) :: y

    real(dp) :: o_a
    integer :: n, incx, incy

    if (present(a)) then
        o_a = a
    else
        o_a = 1.0_dp
    end if

    incx = 1
    incy = 1

    n = size(x)

    call daxpy(n, o_a, x(1), incx, y(1), incy)

end subroutine axpy

!------------------------------------------------------------------------------

subroutine gemm(a, b, c, transa, transb, alpha, beta)

    external :: dgemm

    real(dp), intent(in), optional :: alpha, beta
    character(len=1), intent(in), optional :: transa, transb
    real(dp), dimension(:,:), intent(in) :: a, b
    real(dp), dimension(:,:) , intent(inout) :: c

    integer :: m, n, k, lda, ldb, ldc
    character(len=1) :: o_transa, o_transb
    real(dp) :: o_alpha, o_beta

    if (present(alpha)) then
        o_alpha = alpha
    else
        o_alpha = 1.0_dp
    end if

    if (present(beta)) then
        o_beta = beta
    else
        o_beta = 0.0_dp
    end if

    if (present(transa)) then
        o_transa = transa
    else
        o_transa = 'N'
    end if

    if (present(transb)) then
        o_transb = transb
    else
        o_transb = 'N'
    end if

    if (o_transa == 'N') then
        k = size(a, 2)
    else
        k = size(a, 1)
    end if

    m = size(c, 1)
    n = size(c, 2)
    lda = max(1, size(a, 1))
    ldb = max(1, size(b, 1))
    ldc = max(1, size(c, 1))

    call dgemm(o_transa, o_transb, m, n, k, o_alpha, a(1,1), lda, b(1,1), ldb, o_beta, c(1,1), ldc)

end subroutine gemm

!------------------------------------------------------------------------------

subroutine gemv(a, x, y, alpha, beta, trans)

    external :: dgemv

    real(dp), intent(in), optional :: alpha, beta
    character(len=1), intent(in), optional :: trans
    real(dp), dimension(:,:), intent(in) :: a
    real(dp), dimension(:), intent(in) :: x
    real(dp), dimension(:), intent(inout) :: y

    integer :: m, n, incx, incy, lda
    character(len=1) :: o_trans
    real(dp) :: o_alpha, o_beta

    if (present(alpha)) then
        o_alpha = alpha
    else
        o_alpha = 1.0_dp
    end if

    if (present(beta)) then
        o_beta = beta
    else
        o_beta = 0.0_dp
    end if

    if (present(trans)) then
        o_trans = trans
    else
        o_trans = 'N'
    end if

    incx = 1
    incy = 1
    lda = max(1, size(a, 1))
    m = size(a, 1)
    n = size(a, 2)

    call dgemv(o_trans, m, n, o_alpha, a(1,1), lda, x(1), incx, o_beta, y(1), incy)

end subroutine gemv

!------------------------------------------------------------------------------

subroutine spmv(ap, x, y, uplo, alpha, beta)

    external :: dspmv

    real(dp), dimension(:), intent(in) :: ap
    real(dp), dimension(:), intent(in) :: x
    real(dp), dimension(:), intent(inout) :: y
    character(len=1), intent(in), optional :: uplo
    real(dp), intent(in), optional :: alpha, beta

    integer :: n, incx, incy
    real(dp) :: o_alpha, o_beta
    character(len=1) :: o_uplo

    if (present(uplo)) then
        o_uplo = uplo
    else
        o_uplo = 'U'
    end if

    if (present(alpha)) then
        o_alpha = alpha
    else
        o_alpha = 1.0_dp
    end if

    if (present(beta)) then
        o_beta = beta
    else
        o_beta = 0.0_dp
    end if

    incx = 1
    incy = 1
    n = size(x)

    call dspmv(o_uplo, n, o_alpha, ap(1), x(1), incx, o_beta, y(1), incy)

end subroutine spmv

!------------------------------------------------------------------------------

end module pe_blas_interfaces
