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

!      here we can collect variable type definitions
!      i'm not so sure about the decimals and exponent ranges so please edit
!
!      radovan bast (2008-04-05)
!      inspiration: discussions with jonas juselius
!-------------------------------------------------------------------------------

module kinds_m

  implicit none

  integer, public, parameter :: kind_i2 = selected_int_kind(4)
  integer, public, parameter :: kind_i4 = selected_int_kind(8)

  integer, public, parameter :: kind_sp = selected_real_kind(6, 30)
  integer, public, parameter :: kind_dp = selected_real_kind(14, 200)

end module kinds_m
