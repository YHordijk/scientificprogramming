!      Copyright (c) 2019 by the authors of DIRAC.
!      All Rights Reserved.
!
!      This source code is part of the DIRAC program package.
!      It is provided under a written license and may be used,
!      copied, transmitted, or stored only in accordance to the
!      conditions of that written license.
!
!      In particular, no part of the source code or compiled modules may
!      be distributed outside the research group of the license holder.
!      This means also that persons (e.g. post-docs) leaving the research
!      group of the license holder may not take any part of Dirac,
!      including modified files, with him/her, unless that person has
!      obtained his/her own license.
!
!      For information on how to get a license, as well as the
!      author list and the complete list of contributors to the
!      DIRAC program, see: http://www.diracprogram.org

! Contains a module for elementary 3-vector operations.
! 
! JHS, 17-04-2008
!      

module vector_fun
  implicit none      

  private
  public cross,norm

contains
  function cross(vector1,vector2)
    !     Returns the crossproduct between two vectors
    real(kind=8) :: cross(3)
    real(kind=8) :: vector1(3),vector2(3)
    cross(1)=vector1(2)*vector2(3)-vector1(3)*vector2(2)
    cross(2)=vector1(3)*vector2(1)-vector1(1)*vector2(3)
    cross(3)=vector1(1)*vector2(2)-vector1(2)*vector2(1)
  end function

  real(kind=8) function norm(vector1)
    !     Returns the norm of a vector
    real(kind=8) :: vector1(:)
    norm=sqrt(dot_product(vector1,vector1))
  end function norm

end module
