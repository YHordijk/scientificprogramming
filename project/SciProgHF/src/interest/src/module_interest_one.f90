! ------------------------------------------------------------------------------------
!
! Program:      Dirac 
!
! Library:      InteRest
!
! Module:       module_interest_one.f90 
!
! Description:  InteRest library routine 
!
! Contains:     
!
!                     Copyright (c) 2016 by the authors of DIRAC.
!                     All Rights Reserved.
!                     
!                     This source code is part of the DIRAC program package.
!                     It is provided under a written license and may be used,
!                     copied, transmitted, or stored only in accordance to the
!                     conditions of that written license.
!                     
!                     In particular, no part of the source code or compiled modules may
!                     be distributed outside the research group of the license holder.
!                     This means also that persons (e.g. post-docs) leaving the research
!                     group of the license holder may not take any part of Dirac,
!                     including modified files, with him/her, unless that person has
!                     obtained his/her own license.
!                     
!                     For information on how to get a license, as well as the
!                     author list and the complete list of contributors to the
!                     DIRAC program, see: http://dirac.chem.vu.nl
!
! Author:       Michal Repisky (michal.repisky@uit.no)
!
! Revision:     2.0   
!
! ------------------------------------------------------------------------------------
  subroutine interest_initialize( info )

    use module_interest_osr
    implicit none
    !> input
    logical :: info

    if( is_interest_initialized )then
      return
     !write(6,*)
     !write(6,'(2x,a)')'Integral module InteRest was already initialized'
     !write(6,*)
    else
      if( info )then
        write(6,*)
        write(6,'(2x,a)')'Integral module InteRest is initialized'
        call print_interest_git_revision_info
        write(6,*)
      endif
      call interest_osr_initialize()
      is_interest_initialized = .true.
    endif

  end subroutine
