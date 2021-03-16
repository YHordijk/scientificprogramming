!==============================================================================!
subroutine lap_driver(laplace)
!------------------------------------------------------------------------------!
!
! driver routine that
!
!  computes the numerical Laplace transformation of the orbital energy 
!  denominator by using the Remez algorithm to compute the minimax solution.
!
! the external laplace-minimax routine is called that is hosted on GitHub
!
! BHP, summer 2016
!
!------------------------------------------------------------------------------!

 use type_laplace, only : laplace_info, allocate_laplace

 implicit none

! parameter:
#include "lapdim.h"

! explicit interface
#include "laplace_minimax.h"

! common block:
#include "dcblap.h"

! constants:
 logical, parameter :: locdbg = .false.

!output:
 type(laplace_info), intent(inout) :: laplace

! local:
 real(8) :: errmax

 if (locdbg) write(6,"(a)") "entered lap_driver .."

 call qenter("lap_driver")

 ! determine Laplace quadrature parameters with the minimax algorith
 call laplace_minimax(errmax,xpnts,wghts,nlap,&
                      laplace%bounds(1),laplace%bounds(2),&
                      mxiter=mxiter,iprint=iprlap,stepmx=stepmx,tolrng=tolrng,&
                      tolpar=tolpar,tolerr=tollap,do_rmsd=.false.,&
                      do_init=inilap,do_nlap=fndlap)

 ! allocate Laplace type
 call allocate_laplace(laplace,nlap)

 laplace%num       = nlap
 laplace%errmax    = errmax

 call dcopy(nlap,xpnts,1,laplace%val,1)
 call dcopy(nlap,wghts,1,laplace%wei,1)

 call qexit("lap_driver")

 if (locdbg) write(6,"(a)") "left lap_driver .."

end subroutine lap_driver
!==============================================================================!
