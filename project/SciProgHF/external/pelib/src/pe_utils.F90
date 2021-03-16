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
module pe_utils

    use pe_precision

    implicit none

contains

!------------------------------------------------------------------------------

function elem2charge(elem) result(charge)

    use pe_variables, only: luout

    character(len=*), intent(in) :: elem

    integer :: i
    real(dp) :: charge
    character(len=2), dimension(109) :: elements

    elements = (/ 'H ', 'He', 'Li', 'Be', 'B ', 'C ', 'N ', 'O ', 'F ', 'Ne',&
                & 'Na', 'Mg', 'Al', 'Si', 'P ', 'S ', 'Cl', 'Ar', 'K ', 'Ca',&
                & 'Sc', 'Ti', 'V ', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',&
                & 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y ', 'Zr',&
                & 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn',&
                & 'Sb', 'Te', 'I ', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd',&
                & 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb',&
                & 'Lu', 'Hf', 'Ta', 'W ', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg',&
                & 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th',&
                & 'Pa', 'U ', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm',&
                & 'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt' /)

    if (elem == 'X') then
        charge = 0.0_dp
        return
    end if

    do i = 1, 112
        if (elem == trim(elements(i))) then
            charge = real(i, dp)
            exit
        else
            charge = 0.0_dp
        end if
    end do

    if (int(charge) == 0 .and. elem /= 'X') then
        write(luout, *) 'WARNING: charge not found for element: ', elem
    end if

end function elem2charge

!------------------------------------------------------------------------------

function charge2elem(charge) result(elem)

    real(dp), intent(in) :: charge

    character(len=2) :: elem
    character(len=2), dimension(109) :: elements

    elements = (/ 'H ', 'He', 'Li', 'Be', 'B ', 'C ', 'N ', 'O ', 'F ', 'Ne',&
                & 'Na', 'Mg', 'Al', 'Si', 'P ', 'S ', 'Cl', 'Ar', 'K ', 'Ca',&
                & 'Sc', 'Ti', 'V ', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',&
                & 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y ', 'Zr',&
                & 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn',&
                & 'Sb', 'Te', 'I ', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd',&
                & 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb',&
                & 'Lu', 'Hf', 'Ta', 'W ', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg',&
                & 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th',&
                & 'Pa', 'U ', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm',&
                & 'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt' /)

    if ((nint(charge) == 0) .or. (nint(charge) > 109)) then
        elem = 'X '
        return
    end if

    elem = elements(nint(charge))

end function charge2elem

!------------------------------------------------------------------------------

function charge2mass(charge) result(mass)

    use pe_variables, only: luout

    real(dp), intent(in) :: charge

    integer :: i
    real(dp) :: mass
    real(dp), dimension(109) :: masses

    masses = (/   1.007825_dp,   4.002603_dp,   7.016005_dp,   9.012183_dp,&
              &  11.009305_dp,  12.000000_dp,  14.003074_dp,  15.994915_dp,&
              &  18.998403_dp,  19.992439_dp,  22.989770_dp,  23.985045_dp,&
              &  26.981541_dp,  27.976928_dp,  30.973763_dp,  31.972072_dp,&
              &  34.968853_dp,  39.962383_dp,  38.963708_dp,  39.962591_dp,&
              &  44.955914_dp,  47.947947_dp,  50.943963_dp,  51.940510_dp,&
              &  54.938046_dp,  55.934939_dp,  58.933198_dp,  57.935347_dp,&
              &  62.929599_dp,  63.929145_dp,  68.925581_dp,  73.921179_dp,&
              &  74.921596_dp,  79.916521_dp,  78.918336_dp,  83.911506_dp,&
              &  84.911800_dp,  87.905625_dp,  88.905856_dp,  89.904708_dp,&
              &  92.906378_dp,  97.905405_dp,  97.907216_dp, 101.904348_dp,&
              & 102.905503_dp, 105.903475_dp, 106.905095_dp, 113.903361_dp,&
              & 114.903875_dp, 119.902199_dp, 120.903824_dp, 129.906229_dp,&
              & 126.904477_dp, 131.904148_dp, 132.905429_dp, 137.905232_dp,&
              & 138.906347_dp, 139.905433_dp, 140.907647_dp, 141.907719_dp,&
              & 144.912743_dp, 151.919728_dp, 152.921225_dp, 157.924019_dp,&
              & 158.925342_dp, 163.929171_dp, 164.930319_dp, 165.930290_dp,&
              & 168.934212_dp, 173.938859_dp, 174.940770_dp, 179.946546_dp,&
              & 180.947992_dp, 183.950928_dp, 186.955744_dp, 191.961467_dp,&
              & 192.962917_dp, 194.964766_dp, 196.966543_dp, 201.970617_dp,&
              & 204.974401_dp, 207.976627_dp, 208.980374_dp, 208.982404_dp,&
              & 209.987126_dp, 222.017571_dp, 223.019733_dp, 226.025403_dp,&
              & 227.027750_dp, 232.038051_dp, 231.035880_dp, 238.050785_dp,&
              & 237.048168_dp, 244.064199_dp, 243.061373_dp, 247.070347_dp,&
              & 247.070300_dp, 251.079580_dp, 252.082944_dp, 257.095099_dp,&
              & 258.098570_dp, 259.100931_dp, 260.105320_dp, 261.108690_dp,&
              & 262.113760_dp, 263.118220_dp, 262.122930_dp, 269.134100_dp,&
              & 267.138000_dp /)

    if (int(charge) == 0) then
        mass = 0.0_dp
        return
    end if

    do i = 1, 112
        if (int(charge) == i) then
            mass = masses(i)
            exit
        else
            mass = 0.0_dp
        end if
    end do

    if (int(mass) == 0 .and. int(charge) == 0) then
        write(luout, *) 'WARNING: mass not found for element with charge: ', charge
    end if

end function charge2mass

!------------------------------------------------------------------------------

subroutine chcase(string, uplo)

    character(len=*), intent(inout) :: string
    character(len=*), intent(in), optional :: uplo

    integer :: i, gap
    character(len=1) :: a, z, o_uplo

    if (present(uplo)) then
        o_uplo = uplo
    else
        o_uplo = 'u'
    end if

    gap = iachar('a') - iachar('A')

    if (o_uplo == 'u' .or. o_uplo == 'U') then
        a = 'a'
        z = 'z'
    else if (o_uplo == 'l' .or. o_uplo == 'L') then
        a = 'A'
        z = 'Z'
        gap = - gap
    else
        stop 'Unknown case specified'
    end if

    do i = 1, len_trim(string)
        if (lge(string(i:i), a) .and. lle(string(i:i), z)) then
            string(i:i) = achar(iachar(string(i:i)) - gap)
        end if
    end do

end subroutine chcase

!------------------------------------------------------------------------------

subroutine openfile(filename, lunit, stat, frmt)

    character(*), intent(in) :: filename, stat, frmt
    integer, intent(out) :: lunit
    integer :: i
    logical :: lexist, lopen

    if (stat == 'old') then
        inquire(file=filename, exist=lexist)

        if (.not. lexist) then
            print *, filename, ' not found.'
            stop 'ERROR: file not found.'
        end if
    end if

    do i = 21, 99
        inquire(unit=i, opened=lopen)
        if (lopen) then
            cycle
        else
            lunit = i
            open(unit=lunit, file=filename, status=stat, form=frmt)
            exit
        end if
    end do

    return

end subroutine openfile

!------------------------------------------------------------------------------

end module pe_utils
