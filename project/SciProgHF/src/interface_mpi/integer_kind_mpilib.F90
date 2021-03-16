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
! purpose: store default integer size in Dirac and MPI library
!
module integer_kind_mpilib

#ifdef VAR_MPI

#if defined (INT_STAR8)
#define my_MPI_INTEGER MPI_INTEGER8
#else
#define my_MPI_INTEGER MPI_INTEGER4
#endif

#endif

#ifdef INT_STAR8
  integer, parameter, public :: integer_kind        = 8
#else
  integer, parameter, public :: integer_kind        = 4
#endif
  integer, public            :: integer_kind_in_mpi = 4

!
! next variable is here so it is available both in f77 and f90 interface routines /hjaaj
!
  logical, public            :: MPI_INIT_called     = .false.

end module integer_kind_mpilib
