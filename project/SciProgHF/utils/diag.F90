!************************************************************************ 
!
! Standalone program for comprehensive testing of all diagonalization
! routines used by DIRAC. 
!
! These are Householder, RSJACOBI, QJACOBI and
! DSYEVR.
!
! The program either reads the matrix file, or generates a random matrix.
!
!
! Written by Miro Ilias, 07,08-2014 Banska Bystrica, Slovakia
! with the help of P.S. Bagus and T.Saue
!
!************************************************************************ 

#ifdef MOD_UNRELEASED
module time
  implicit none
  character*10 start_date, end_date ! format of date: "dd_mmm_yy "
  character*8  start_time, end_time ! format of time: "hh:mm:ss"
end module time
#endif /* MOD_UNRELEASED */

Program Test_DIRAC_Diagonalization_Routines
#ifdef MOD_UNRELEASED
 use diagmod
#include "priunit.h"
  write(LUPRI,'(/,2x,a,/)') "**** Welcome to the testing program of DIRAC diagonalization routines ! ****"
  call DBLGRP !> initialize include file variable for printing matrices (PRQMAT)
  call read_input_file
  call read_matrix
  if (do_rs)      call test_dirac_rs
  if (do_dsyevr)  call test_lapack_dsyevr
  if (do_dsyevr_ulf)  call test_lapack_dsyevr_ulf
  if (do_qjaco)   call QJACOBITEST
  if (do_householder_psb) call test_PSB_householder
#else
 call quit('program not in release version!')
#endif /* MOD_UNRELEASED */
End Program Test_DIRAC_Diagonalization_Routines

!--------------------------------------------------------------
! placed this routine here as it depends on src/input module
! can not be in src/gp module !
!--------------------------------------------------------------
subroutine read_input_file
  implicit none
! reads the text input file "DIAGONALIZE_TESTS.INP"
  character*21 :: inpfile="DIAGONALIZE_TESTS.INP"
  logical :: input_found
#include "priunit.h"
    inquire(file=inpfile, exist=input_found)
    if (.not.input_found) then
      write(LUPRI,"(2x,a,a)") "read_input_file: the input file not found ! input_file=",inpfile
      call quit("read_input_file: the input file not found !")
    endif
    call read_menu_input(inpfile, LUSTDIN, 'ALL', input_found)
end subroutine read_input_file

