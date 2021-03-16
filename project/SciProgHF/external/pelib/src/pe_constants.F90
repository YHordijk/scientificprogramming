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
module pe_constants

    use pe_precision

    ! 1 bohr = 0.5291772109217 Aa (codata 2010)
    real(dp), parameter :: bohr2aa = 0.5291772109217_dp
    real(dp), parameter :: aa2bohr = 1.0_dp / bohr2aa
    real(dp), parameter :: kcal2hartree = 1.0_dp / 627.5091809_dp
    real(dp), parameter :: pi = 3.141592653589793_dp
    real(dp), parameter :: hartree2ev = 27.211_dp
    real(dp), parameter :: ev2cmi = 8065.73_dp
    real(dp), parameter :: hartree2cmi = hartree2ev * ev2cmi

end module pe_constants
