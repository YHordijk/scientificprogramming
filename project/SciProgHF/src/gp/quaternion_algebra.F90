module quaternion_algebra

!  This module provides multiplication rules for
!  quaternion unit vectors 1,i,j and k
!
!  IQMULT(e_1,e_2,e_M) and IQPHAS(e_1,e_2,e_M) gives the
!  resulting unit vector and phase (+/- 1) of a quaternion
!  triple product (e_1)(e_M)(e_2) :
!
!  0 |  1  i  j  k    i |  1  i  j  k
!  ---------------    ---------------
!  1 |  1  i  j  k    1 |  i -1  k -j
!  i |  i -1  k -j    i | -1 -i -j -k
!  j |  j -k -1  i    j | -k -j  i  1
!  k |  k  j -i -1    k |  j -k -1  i
!
!  j |  1  i  j  k    k |  1  i  j  k
!  ---------------    ---------------
!  1 |  j -k -1  i    1 |  k  j -i -1
!  i |  k  j -i -1    i | -j  k  1 -i
!  j | -1 -i -j -k    j |  i -1  k -j
!  k | -i  1 -k  j    k | -1 -i -j -k
!
!  IQSIGN(e_1,FX,TX) gives the sign of the quaternion units under
!                    various operations:
!  FX = 1 : Standard quaternion:   q  = a  + bj
!  FX = 2 : Complex conjugation:   q* = a* - bj
!
!  TX = 1 : Standard quaternion:   q  = a  + bj
!  TX = 2 : i-transform        : -iqi = a  - bj
!  TX = 3 : j-transform        : -jqj = a* + b*j
!  TX = 4 : k-transform        : -kqk = a* - b*j

   implicit none

   integer, parameter :: iqmult(4, 4, 4) = reshape((/1, 2, 3, 4, 2, 1, 4, 3, 3, 4, 1, 2, 4, 3, 2, 1,   &
                                                     2, 1, 4, 3, 1, 2, 3, 4, 4, 3, 2, 1, 3, 4, 1, 2,   &
                                                     3, 4, 1, 2, 4, 3, 2, 1, 1, 2, 3, 4, 2, 1, 4, 3,   &
                                                     4, 3, 2, 1, 3, 4, 1, 2, 2, 1, 4, 3, 1, 2, 3, 4/), &
                                                   (/4, 4, 4/))

   integer, parameter :: IQPHASE(4, 4, 4) = reshape((/1, 1, 1, 1, 1,-1,-1, 1, 1, 1,-1,-1, 1,-1, 1,-1,   &
                                                      1,-1,-1, 1,-1,-1,-1,-1, 1,-1, 1,-1,-1,-1, 1, 1,   &
                                                      1, 1,-1,-1,-1, 1,-1, 1,-1,-1,-1,-1, 1,-1,-1, 1,   &
                                                      1,-1, 1,-1, 1, 1,-1,-1,-1, 1, 1,-1,-1,-1,-1,-1/), &
                                                    (/4, 4, 4/))

   integer, parameter :: iqsign(4, 2, 4) = reshape((/1, 1, 1, 1, 1,-1,-1,-1,   &
                                                     1, 1,-1,-1, 1,-1, 1, 1,   &
                                                     1,-1, 1,-1, 1, 1,-1, 1,   &
                                                     1,-1,-1, 1, 1, 1, 1,-1/), &
                                                   (/4, 2, 4/))

   character(1), parameter :: qunit(4) = (/'1','i','j','k'/)

end module
