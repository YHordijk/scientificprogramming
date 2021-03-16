!-----------------------------------------------------------------------------------------
!
! Collection of Fortran90 modules encapsulating the "old" Fortran77 include files
! (dcbham.h,dcbgen.h,...) in the Dirac source codes.
!
!  All variables inside the h-include files become PRIVATE within the dedicated Fotran90 module and are
! accessible out of module scope ONLY via dedicated subroutines/functions (methods).
!
! Programmers can still use the h-include files in the old way (#include "dcbgham.h")
!
!
! Notation examples:  
! ------------------
!
!             use include_dcbgen_h 
!             use include_dcbham_h
!             use include_dgroup_h, only : get_IPQTOQ
!             use include_dcbdhf_h
!
! Examples of usage:
! ------------------
!
!    Function get_<variable>(parameter1,parameter2,...) - function returns include file (module) variable
!
!    Subroutine set_<variable>(parameter1,parameter2,...) - method changes include file (module) variable
!
!
! TODO: extend encapsulating modules for other inlcude files and other variables
!
!
! Written by Miro Ilias, August 2014, Banska Bystrica; July 2016, GSI.de, Darmstadt
!
!-----------------------------------------------------------------------------------------

module include_dcbgen_h
implicit none
private
#include "dcbgen.h"
!> specify public variables and public methods
public :: get_IPRGEN
contains
function get_IPRGEN() result (IPRGEN_out)
 implicit none
 integer :: IPRGEN_out
 IPRGEN_out = IPRGEN
end function get_IPRGEN
end module include_dcbgen_h

module include_dgroup_h
!> module for picking up/setting certain variables from dgroup.h
implicit none
private
#include "dgroup.h"
!> specify public variables and public methods to operate upon internal (protected) data module
public :: get_IPQTOQ
public :: get_NZ  
contains
function get_IPQTOQ(i,j) result (IPQTOQ_out)
!> get protected (internal) indexed variable IPQTOQ(i,j) from the include file "dgroup.h"
 implicit none
 integer, intent(in) :: i,j
 integer :: IPQTOQ_out
 if (.not. (i>=1.and.i<=4.and.j>=0.and.j<=7)) then
   call quit('get_IPQTOQ(i,j): wrong values i,j !')
 endif
 IPQTOQ_out=IPQTOQ(i,j)
end function get_IPQTOQ

function get_NZ() result (NZ_out)
!> retrieves the important NZ value from "dgroup.h"
 implicit none
 integer :: NZ_out
 NZ_out=NZ
end function get_NZ

end module include_dgroup_h

module include_dcbham_h
!> module encapsulating the dcbham.h include file
!> can pick up/set up some defined variables in dcbham.h
implicit none
private
#include"dcbham.h"
!> specify public variables and public methods
public :: get_IPRHAM
contains
function get_IPRHAM() result (IPRHAM_out)
 implicit none
 integer :: IPRHAM_out
 IPRHAM_out = IPRHAM
end function get_IPRHAM
end module include_dcbham_h

module include_dcbdhf_h
!> module encapsulating the "dcbdhf.h" include file
implicit none
private
#include "dcbdhf.h"
!> specify public variables and public methods
public :: get_IPRSCF
contains
function get_IPRSCF() result (IPRSCF_out)
 implicit none
 integer :: IPRSCF_out
 IPRSCF_out = IPRSCF
end function get_IPRSCF
end module include_dcbdhf_h

module include_dcbprp_h
!> module encapsulating the "dcbprp.h" include file
implicit none
private
!> needed for parameters used in dcbprp.h
#include "mxcent.h"
!> own private include file
#include "dcbprp.h"
!> specify public variables and public methods
public :: get_IPRPRP
contains
function get_IPRPRP() result (IPRPRP_out)
 implicit none
 integer :: IPRPRP_out
 IPRPRP_out = IPRPRP
end function get_IPRPRP
end module include_dcbprp_h
