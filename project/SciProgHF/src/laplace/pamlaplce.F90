!==============================================================================!
subroutine pamlaplce()
!------------------------------------------------------------------------------!
!
!
!                              L A P L A C E
!
!
!                uses the numerical Laplace transformation for
!                       orbital energy denominators
!
!
!  author:         Benjamin Helmich-Paris
!
!  variables taken from EVIL EVIL common blocks:
!
!------------------------------------------------------------------------------!
 
 use type_laplace, only : laplace_info, deallocate_laplace

 implicit none

! parameters
#include "lapdim.h"

! common blocks:
#include "dcblap.h"
#include "dcbgen.h"
#include "dcborb.h"

! local parameter:
 character(len=*),parameter :: chrdbg     = "pamlaplce>"

! local scalars:
 type(laplace_info) :: laplace

 integer :: noct, nvirct, nfro_lap, nfrv_lap
 integer :: idum
 integer :: istro,istrv,iendo, iendv
 logical :: filex
 real(8) :: dum, toterg

! local allocatables:
 real(8), allocatable :: eig(:)

 call qenter("pamlaplce")

!------------------------------------------------------------------------------!
! echo AOMP2 input parameter
!------------------------------------------------------------------------------!
 write(6,'(/3a/)') '+',repeat('=',70),'+'
 write(6,'(/2a/)')     repeat(' ',27),'L A P L A C E   module'
 write(6,'(/3a/)') '+',repeat('=',70),'+'


!------------------------------------------------------------------------------!
! orbital occupation 
!
! @trond: Please figure out the number of frozen occupied and virtual orbitals
!------------------------------------------------------------------------------!
 nfro_lap = 0
 nfrv_lap = 0

 noct   = nish(1)-nfro_lap
 nvirct = nssh(1)-nfrv_lap

 write(6,'(////,a,/)') "         orbital occupation                "
 write(6,'(a,/)')      "  -------------------------------          "
 write(6,'(a,i5,/)')   "   positronic MOs: ",npsh(1)
 write(6,'(a,i5,/)')   "   electronic MOs: ",nesh(1)
 write(6,'(a,i5)')     "   occupied: ",noct+nfro_lap
 if (nfro_lap.gt.0) then
  write(6,'(a,i5)')    "    frz. occupied: ",nfro_lap
  write(6,'(a,i5,/)')  "    act. occupied: ",noct
 end if
 write(6,'(a,i5)')     "   virtual: ",nvirct+nfrv_lap
 if (nfrv_lap.gt.0) then
  write(6,'(a,i5)')    "    act.  virtual: ",nvirct
  write(6,'(a,i5)')    "    frz.  virtual: ",nfrv_lap
 end if
 write(6,'(a)')        "  -------------------------------          "
 write(6,'(a,i5,//)')  "   total # of MOs: ",norb(1)
 !write(6,'(a,i5,//)')  "   total # of cart. AOs: ",nbfcao
 !write(6,'(a,i5,//)')  "   total # of sphr. AOs: ",nbfsao

!------------------------------------------------------------------------------!
! Laplace quadrature 
!------------------------------------------------------------------------------!

  allocate(eig(npsh(1)+nesh(1)))

  ! read MO energies from DFCOEF
  INQUIRE(FILE='DFCOEF',EXIST=FILEX)
  IF(.NOT.FILEX) CALL QUIT ('DFCOEF missing')
  CALL OPNFIL(LUCOEF,'DFCOEF','OLD','PRPEXP')
  call reacmo(LUCOEF, 'DFCOEF', dum, eig, idum, toterg, 4  )
  CLOSE(LUCOEF,STATUS='KEEP')

  ! orbital start and end addresses
  istro = npsh(1)+nfro_lap+1
  iendo = npsh(1)+nfro_lap+noct
  istrv = npsh(1)+nfro_lap+noct+1
  iendv = npsh(1)+nfro_lap+noct+nvirct
 
  ! bounds of orbital energy denominator
  laplace%bounds(1) = 2.d0*(eig(istrv)-eig(iendo))
  laplace%bounds(2) = 2.d0*(eig(iendv)-eig(istro))

  deallocate(eig)

  call lap_driver(laplace)

 call deallocate_laplace(laplace)

 call qexit("pamlaplce")
 
end subroutine pamlaplce
!==============================================================================!
