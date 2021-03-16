  subroutine localization_interface
  !==============================================================================!
  !           Interface of the minimization module with Dirac program            !
  !==============================================================================!
  implicit none

#include "dcbloc.h"
#include "priunit.h"
#include "dcbgen.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dcbham.h"
#include "dgroup.h"
#include "mxcent.h"
#include "nuclei.h"
#include "thrzer.h"
#include "dcbprj.h"

  ! Internal variables
  integer :: i, ier, iunit, ifrp

  integer, allocatable :: indao(:), indaos(:), indnao(:)
  real*8,  allocatable :: smat(:)

  !------------------------------------------------------------------------------!
  !-----------------------------------------------------!
  ! Interface with the second_order_minimization module !
  !-----------------------------------------------------!
  iunit = 60

  ! open file
  open(iunit,file='OPTIM.INP',status='new',form='formatted',access='sequential',iostat=ier)

  if(ier /= 0) call quit('localization_interface: ERROR while creating new OPTIM.INP file')
  rewind iunit

  write(iunit,*) 'iprloc  ',iprloc
  write(iunit,*) 'itrloc  ',itrloc
  write(iunit,*) 'lupri   ',lupri
  write(iunit,*) 'thfull  ',thfull
  write(iunit,*) 'thdiag  ',thdiag
  write(iunit,*) 'thgrad  ',thgrad
  write(iunit,*) 'lgfull  ',lgfull
  write(iunit,*) 'lgdiag  ',lgdiag
  write(iunit,*) 'lggrad  ',lggrad
  write(iunit,*) 'lgchck  ',lgchck
  write(iunit,*) 'hesloc  ',hesloc

  ! close file
  close(iunit,status='keep')

  !-------------------------------------------!
  ! Interface with the func_PipekMezey module !
  !-------------------------------------------!
  ! Calculate indnao
  ifrp = 1
  allocate( indao (nfbas(ifrp,0)), stat=ier )
  if(ier /= 0)then
    write(lupri,*) 'ERROR: Problem with allocation of indao in localization_interface'
    stop
  endif
  allocate( indaos(nfbas(ifrp,0)), stat=ier )
  if(ier /= 0)then
    write(lupri,*) 'ERROR: Problem with allocation of indaos in localization_interface'
    stop
  endif
  allocate( indnao(2*nucdep+2),    stat=ier )
  if(ier /= 0)then
    write(lupri,*) 'ERROR: Problem with allocation of indnao in localization_interface'
    stop
  endif
  call indaoc(indao,indaos,indnao)

  ! open file
  open(iunit,file='FUNC_PM.INP',status='new',form='formatted',access='sequential',iostat=ier)

  if(ier /= 0) call quit('localization_interface: ERROR while creating new FUNC_PM.INP file')
  rewind iunit

  ! dcbloc.h, priunit.h
  write(iunit,*) 'iprloc  ',iprloc
  write(iunit,*) 'prjloc  ',prjloc
  write(iunit,*) 'lupri   ',lupri
  write(iunit,'(A72)') selmos
  ! dcbgen.h
  write(iunit,*) title
  write(iunit,*) 'lucoef  ',lucoef
  ! dcbbas.h
  write(iunit,*) 'n2bbasx ',n2bbasx
  write(iunit,*) 'n2bbasxq',n2bbasxq
  write(iunit,*) 'nfbas   '
  write(iunit,*) nfbas(:,:)
  ! dcborb.h
  write(iunit,*) 'i_dcborb',i_dcborb_set
  write(iunit,*) 'norbt   ',norbt
  write(iunit,*) 'ncmotq  ',ncmotq
  write(iunit,*) 'nish    ',nish(:)
  write(iunit,*) 'norb    ',norb(:)
  write(iunit,*) 'npsh    ',npsh(:)
  write(iunit,*) 'nesh    ',nesh(:)
  write(iunit,*) 'iorb    ',iorb(:)
  write(iunit,*) 'icmoq   ',icmoq(:)
  ! dcbham.h
  write(iunit,*) 'ssmtrc  ',ssmtrc
  write(iunit,*) 'levyle  ',levyle
  ! dgroup.h
  write(iunit,*) 'nfsym   ',nfsym
  write(iunit,*) 'nz      ',nz
  write(iunit,*) 'ipqtoq  '
  write(iunit,*) ipqtoq(:,:)
  ! nuclei.h
  write(iunit,*) 'nucdep  ',nucdep
  ! thrzer.h
  write(iunit,*) 'thrzer  ',thrzer
  ! dcbprj.h
  write(iunit,*) 'nrefs   ',nrefs
  write(iunit,*) 'iprprj  ',iprprj
  write(iunit,*) 'prothr  ',prothr
  write(iunit,*) 'ownbas  ',ownbas
  write(iunit,*) 'maxref  ',maxref
  write(iunit,*) 'twocomp ',twocomp
  write(iunit,*) 'x2c     ',x2c
  do i=1,nrefs
    write(iunit,*) reffil(i)
    write(iunit,*) vecref(1,i)
  enddo
  ! indaoc routine
  write(iunit,*) 'indnao  '
  write(iunit,*) indnao(:)

  if(nfsym /= 1)then
    write(lupri,*) 'ERROR: NFSYM /= 0; localization work only with C1 symmetry'
    stop
  endif

  ! close file
  close(iunit,status='keep')

  !----------------------------------!
  ! Write overlap matrix to the file !
  !----------------------------------!
  allocate( smat(n2bbasx), stat=ier )
  if(ier /= 0)then
    write(lupri,*) 'ERROR: Problem with allocation of smat in localization_interface'
    stop
  endif

  call gtovlx(smat,ssmtrc)

  ! open file
  open(iunit,file='OVERLAP.INP',status='new',form='formatted',access='sequential',iostat=ier)

  if(ier /= 0) call quit('localization_interface: ERROR while creating new OVERLAP.INP file')
  rewind iunit

  write(iunit,*) smat

  ! close file
  close(iunit,status='keep')

  deallocate( smat   )
  deallocate( indao  )
  deallocate( indaos )
  deallocate( indnao )

  end subroutine localization_interface
