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
!
! purpose: determine integer size in MPI library to 
!
module integer_model_in_mpi
#ifdef VAR_MPI
  use integer_kind_mpilib
#ifdef USE_MPI_MOD_F90
  use mpi
  implicit none
#else
  implicit none
#include "mpif.h"
#endif

  public determine_integer_model_in_mpi
  public get_integer_model_in_mpi
  
  private
  save

contains

  subroutine determine_integer_model_in_mpi

#ifdef VAR_MPI_32BIT_INT
  integer_kind_in_mpi = 4
#else
  integer_kind_in_mpi = sizeof(MPI_INTEGER)
#endif

  end subroutine determine_integer_model_in_mpi

  function get_integer_model_in_mpi() result(int_mpi)
    integer :: int_mpi
    int_mpi = integer_kind_in_mpi
  end function get_integer_model_in_mpi

#else
  implicit none
  integer         :: blabla_integer_model_in_mpi
#endif
end module integer_model_in_mpi
