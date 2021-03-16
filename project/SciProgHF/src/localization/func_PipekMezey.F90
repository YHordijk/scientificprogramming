module func_PipekMezey
!================================================================================!
!!! This is not a true module, because it still depends on external subroutines!!!
!                                                                                !
!                  This module is used in optimization module                    !
!           it will provide information about Pipek-Mezey functional             !
!                                                                                !
!   Two forms of the Pipek-Mezey functional are available:                       !
!   1) based on Mulliken population analysis  -->  orblo1                        !
!   2) based on projection analysis           -->  ORBLP1                        !
!                                                                                !
!   [1] S. Dubillard, J.-B. Rota, and T. Saue, J.Chem.Phys. 124, 154307 (2006)   !
!================================================================================!
! Input variables of the func_PipekMezey module:                                 !
!                                                                                !
!   prjloc - localization based on:                                              !
!            True  - projection analysis                                         !
!            False - Mulliken analysis                                           !
!   iprloc - print level (no prints =0)                                          !
!   isel   - integer array which contains numbers of MOs that will be localized  !
!   ... ( not finished yet )                                                     !
!================================================================================!
implicit none

public  :: init_func, quit_func, mv2point, get_func, get_gradient, get_hessian, get_date

private :: ORBLP1, ORBLP2, orblo1, ORBLO4, wcmo2file, wricmo_mod, reacmo_mod
private :: ROTCMO_MOD, NUMLST_MOD, errall, get_overlap, selcfs_mod, qtrans90_mod, qdiag90_mod

type :: PipekMezey_flags
  integer :: iprloc
  logical :: prjloc
end type PipekMezey_flags

type :: PipekMezey_variales
  integer :: lupri, lucoef
  integer :: nfbas(2,0:2), n2bbasx, n2bbasxq
  integer :: nish(2), norb(2), npsh(2), nesh(2), iorb(2), icmoq(2), norbt, ncmotq, i_dcborb_set
  real*8  :: ssmtrc
  logical :: h2comp, levyle
  integer :: nfsym, nz, ipqtoq(4,0:7)
  integer :: nucdep
  real*8  :: thrzer
  integer :: nrefs, iprprj
  integer :: maxref
  real*8  :: prothr
  logical :: ownbas
  character*50 :: title
  real*8  :: timstr
end type PipekMezey_variales

type( PipekMezey_flags     ), private :: PM
type( PipekMezey_variales  ), private :: VAR

! Some global parameters
real*8, parameter, private :: d0 = 0.0d0, d1 = 1.0d0, d2 = 2.0d0, d4 = 4.0d0

! functional value of our functional (stored for output)
integer, save :: keyform
real*8,  save :: dinvmg4out

! Global arrays
real*8,  allocatable, save, private :: pmat(:,:,:,:), smat(:)
real*8,  allocatable, save, private :: scoef(:,:,:), storec(:,:,:)
integer, allocatable, save, private :: nfnuc(:), isel(:)

character(len=6 ), allocatable, save, private :: VAR_reffil(:)    ! reffil(maxref)
character(len=72), allocatable, save, private :: VAR_vecref(:,:)  ! vecref(2,maxref)

contains

  subroutine init_func(ndim)
  !==============================================================================!
  !            Initialize module func_PipekMezey for Dirac program               !
  !==============================================================================!
  implicit none

  ! External variables
  integer :: ndim

  ! Internal variables

  logical :: tobe, twocomp, x2c
  integer :: i, ifrp, idum, ier, nloc, nmax, iunit, iopt,ioer
  real*8  :: toterg, second
  character*8  :: chdum
  character*72 :: selmos

  integer, allocatable, dimension(:) :: itmp, indnao, ibeig
  real*8,  allocatable, dimension(:) :: cmo, eig

  !------------------------------------------------------------------------------!
  ! Start localization of Molecular orbitals
  call titler('Molecular orbital localization','*',122)

  !------------------------------------------------------------------------------!
  !                       Initialize global variables                            !
  !------------------------------------------------------------------------------!
  iunit = 61

  ! open file
  open(iunit,file='FUNC_PM.INP',status='old',form='formatted',access='sequential',iostat=ier)

  if(ier /= 0)then
    write(6,*) 'init_func: ERROR while opening FUNC_PM.INP file'
    stop
  endif

  rewind iunit

  ! dcbloc.h, priunit.h
  read(iunit,*) chdum,PM%iprloc
  read(iunit,*) chdum,PM%prjloc
  read(iunit,*) chdum,var%lupri
  read(iunit,'(A72)') selmos
  ! dcbgen.h
  read(iunit,*) var%title
  read(iunit,*) chdum,var%lucoef
  ! dcbbas.h
  read(iunit,*) chdum,var%n2bbasx
  read(iunit,*) chdum,var%n2bbasxq
  read(iunit,*) chdum
  read(iunit,*) var%nfbas(:,:)
  ! dcborb.h
  read(iunit,*) chdum,var%i_dcborb_set
  read(iunit,*) chdum,var%norbt
  read(iunit,*) chdum,var%ncmotq
  read(iunit,*) chdum,var%nish(:)
  read(iunit,*) chdum,var%norb(:)
  read(iunit,*) chdum,var%npsh(:)
  read(iunit,*) chdum,var%nesh(:)
  read(iunit,*) chdum,var%iorb(:)
  read(iunit,*) chdum,var%icmoq(:)
  ! dcbham.h
  read(iunit,*) chdum,var%ssmtrc
  read(iunit,*) chdum,var%levyle
  ! dgroup.h
  read(iunit,*) chdum,var%nfsym
  read(iunit,*) chdum,var%nz
  read(iunit,*) chdum
  read(iunit,*) var%ipqtoq(:,:)
  ! nuclei.h
  read(iunit,*) chdum,var%nucdep
  ! thrzer.h
  read(iunit,*) chdum,var%thrzer
  ! dcbprj.h
  read(iunit,*) chdum,var%nrefs
  read(iunit,*) chdum,var%iprprj
  read(iunit,*) chdum,var%prothr
  read(iunit,*) chdum,var%ownbas
  read(iunit,*) chdum,var%maxref
  read(iunit,*) chdum,twocomp
  read(iunit,*) chdum,x2c

  allocate( VAR_reffil(var%maxref),   stat=ier ); call errall(ier,'init_func','VAR_reffil')
  allocate( VAR_vecref(2,var%maxref), stat=ier ); call errall(ier,'init_func','VAR_vecref')

  do i=1,var%nrefs
    read(iunit,*) VAR_reffil(i)
    read(iunit,*) VAR_vecref(1,i)
  enddo

  ! indaoc routine
  allocate( indnao(2*var%nucdep+2), stat=ier ); call errall(ier,'init_func','indnao')
  indnao = 0
  read(iunit,*) chdum
  read(iunit,*) indnao(:)

  ! two component calculation
  var%h2comp = twocomp .or. x2c

  ! starting time
  var%timstr = second()

  ! close file
  close(iunit,status='keep')

  !----------------!
  ! set dimensions !
  !----------------!
  ifrp = 1
  ! set limits on possible localized orbitals
  ! when selmos == 'allocc' subroutine numlst_mod will select 0..nmax MOs
  nmax = var%nesh(ifrp)
  if(selmos == 'allocc') nmax = var%nish(ifrp)
  ! set maximum possible number of localized orbitals
  nloc = var%norb(ifrp)
  ! alocate temporary array; the data will be later moved to isel
  allocate( itmp(nloc), stat=ier ); call errall(ier,'init_func','itmp')
  itmp = 0; idum = 1
  ! read string selmos
  call numlst_mod(selmos,itmp,nloc,0,nmax,ifrp,idum)
  nloc = idum                          ! set number of orbitals for localization
  ndim = 4 * ( nloc*nloc - nloc ) / 2  ! set dimension of your space

  !------------------------------------------------------------------------------!
  !                             Print debug info                                 !
  !------------------------------------------------------------------------------!
  if(PM%iprloc.ge.4)then
    write(var%lupri,*) 'lupri    =',var%lupri
    write(var%lupri,*) 'nish     =',var%nish(ifrp)
    write(var%lupri,*) 'norb     =',var%norb(ifrp)
    write(var%lupri,*) 'nucdep   =',var%nucdep
    write(var%lupri,*) 'maxref   =',var%maxref
    write(var%lupri,*) 'thrzer   =',var%thrzer
    write(var%lupri,*) 'nfbas(0) =',var%nfbas(ifrp,0)
    write(var%lupri,*) 'nfbas(1) =',var%nfbas(ifrp,1)
    write(var%lupri,*) 'nfbas(2) =',var%nfbas(ifrp,2)
    write(var%lupri,*) 'selmos   =',selmos
  endif

  !------------------------------------------------------------------------------!
  !                              Allocate arrays                                 !
  !------------------------------------------------------------------------------!
  ! allocate global arrays
  allocate( scoef (var%nfbas(ifrp,0),nloc,4), stat=ier ); call errall(ier,'init_func','scoef ')
  allocate( storec(var%nfbas(ifrp,0),nloc,4), stat=ier ); call errall(ier,'init_func','storec')
  allocate( pmat  (nloc,nloc,4,var%nucdep),   stat=ier ); call errall(ier,'init_func','pmat  ')
  allocate( smat  (var%n2bbasx),              stat=ier ); call errall(ier,'init_func','smat  ')
  allocate( nfnuc (2*var%nucdep),             stat=ier ); call errall(ier,'init_func','nfnuc ')
  allocate( isel  (nloc),                     stat=ier ); call errall(ier,'init_func','isel  ')
  scoef = d0; storec = d0; pmat = d0; smat = d0; nfnuc = 0

  ! allocate local arrays
  allocate( cmo   (var%n2bbasxq),      stat=ier ); call errall(ier,'init_func','cmo   ')
  allocate( eig   (var%norbt),         stat=ier ); call errall(ier,'init_func','eig   ')
  allocate( ibeig (var%norbt),         stat=ier ); call errall(ier,'init_func','ibeig ')
  cmo = d0

  !------------------------------------------------------------------------------!
  !                      Precalculate some global arrays                         !
  !------------------------------------------------------------------------------!
  !----------------------------------------------------!
  ! -- nfnuc -- Get the number of functions per center !
  !----------------------------------------------------!
  call orblo4(indnao,var%nucdep)

  !-------------------------------!
  ! -- smat -- Get overlap matrix !
  !-------------------------------!
  call get_overlap(smat)
  if(PM%iprloc.ge.6)then
    write(var%lupri,'(/a)') '* orbloc: Overlap matrix' 
    call prqmat(smat,var%nfbas(ifrp,0),var%nfbas(ifrp,0),var%nfbas(ifrp,0),var%nfbas(ifrp,0), &
                1,var%ipqtoq(1,0),var%lupri)
  endif

  !-----------------------------------------------------!
  ! -- scoef -- Electronic orbital coeffients selection !
  !-----------------------------------------------------!
  ! Check if coefficients are on file
  keyform = 0
  inquire(file='DFCOEF',exist=tobe)
  if(tobe) keyform = 1
!  inquire(file='DFPCMO',exist=tobe)
!  if(tobe) keyform = 2
  if(keyform == 0)then
    write(var%lupri,'(a)') 'init_func: Coefficient file not found!'
!    write(var%lupri,'(a)') 'init_func: DFCOEF or DFPCMO file is missing!'
    call quit('init_func: Coefficients not found !')
  endif

  ! Read coefficients
  !  call reacmo_mod(var%lucoef,cmo,eig,ibeig,toterg,1)
  open(var%lucoef,file='DFCOEF',status='old',form='unformatted',access='sequential',iostat=ioer)
  if(ioer /= 0) call quit('REACMO_MOD: ERROR while opening DFCOEF or DFPCMO')
  iopt = 2
  call REACMO(var%lucoef,'DFCOEF',cmo,eig,ibeig,toterg,iopt)
  close(var%lucoef,status='keep')
  
  ! select orbitals
  isel(1:nloc) = itmp(1:nloc)
  call selcfs_mod(cmo,ifrp,scoef,nloc,isel,0,nloc,var%nfbas(ifrp,0),var%norb(ifrp))
  if(PM%iprloc.ge.5)then
    write(var%lupri,'(/a)') '* orblo1: Selected molecular coefficients'
    call prqmat(scoef,var%nfbas(ifrp,0),nloc,var%nfbas(ifrp,0),nloc,4,var%ipqtoq(1,0),var%lupri)
  endif

  !------------------------------------------------------------------------------!
  !                         Deallocate local arrays                              !
  !------------------------------------------------------------------------------!
  deallocate( cmo    )
  deallocate( eig    )
  deallocate( ibeig  )
  deallocate( itmp   )
  deallocate( indnao )

  end subroutine init_func

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine quit_func(hmat)
  !==============================================================================!
  !                  Uninitialize module func_PipekMezey                         !
  !==============================================================================!
  implicit none

  ! External variables
  real*8,  intent(inout), pointer :: hmat(:,:)

  ! Internal variables
  character :: resst*50
  real*8    :: timend, second

  ! Put the localised coefficients in DFCOEF (scoef)
  call wcmo2file

  ! deallocate global arrays
  deallocate( storec )
  deallocate( scoef  )
  deallocate( pmat   )
  deallocate( smat   )
  deallocate( nfnuc  )
  deallocate( isel   )

  deallocate( VAR_reffil )
  deallocate( VAR_vecref )

  ! deallocate pointers
  deallocate( hmat   )

  ! print final functional value
  write(resst,'(a,f21.15)') 'Inverse mean delocalization =',dinvmg4out
  call titler(resst,'*',104)

  ! calculate time spend on localization
  timend = second()
  call timtxt('          --- Time used in localization is :',timend-var%timstr,var%lupri)

  end subroutine quit_func

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine mv2point(restore,rejected,xtmp)
  !==============================================================================!
  !                    Move to the new expansion point                           !
  !                      --> get new MO coefficients                             !
  !                                                                              !
  !   Note:  Deffinition of kappa matrix in rotcmo_mod is different as in        !
  !          calculation of gradient and hessian in get_gradient and get_hessian !
  !          --> Indexes are exchanged and j quaterion is with minus sign        !
  !              Qkappa_lk = kappa_kl - bar{kappa}_kl^* j                        !
  !              It follows from definition of bared unitary matrix:             !
  !              tilde{Psi_i} = U_ij Psi_j + bar{U}_ij Psi_j                     !
  !==============================================================================!
  implicit none

  ! Global variables
  ! real*8,  intent(inout) :: scoef(:,:,:), storec(*)
  ! real*8,  intent(  out) :: pmat(nloc,nloc,4,nucdep)

  ! External variables
  logical, intent(in) :: restore, rejected
  real*8,  intent(in) :: xtmp(:)

  ! Internal variables
  integer :: i, j, ij1, ij2, ij3, ij4, iz, nloc, nloc12, ifrp, ier

  real*8, allocatable :: xkappa(:,:,:)

  !------------------------------------------------------------------------------!
  ifrp = 1
  nloc = size(pmat,1)
  nloc12 = size(xtmp,1) / 4

  ! allocate local variable
  allocate( xkappa(nloc,nloc,4), stat=ier ); call errall(ier,'mv2point','xkappa')

  ! Unpack X to kappa matrix
  xkappa(:,:,:) = d0
  ij1 = 0
  do j = 1,nloc
    do i = j+1,nloc
      ij1 = ij1 + 1
      ij2 = ij1 + nloc12
      ij3 = ij2 + nloc12
      ij4 = ij3 + nloc12
      xkappa(j,i,1) =   xtmp(ij1)
      xkappa(i,j,1) = - xkappa(j,i,1)

      xkappa(j,i,2) =   xtmp(ij3)
      xkappa(i,j,2) =   xkappa(j,i,2)

      xkappa(j,i,3) = - xtmp(ij2)
      xkappa(i,j,3) =   xkappa(j,i,3)

      xkappa(j,i,4) =   xtmp(ij4)
      xkappa(i,j,4) =   xkappa(j,i,4)
    enddo
  enddo

  if(PM%iprloc.ge.2)then
    do iz=1,4
      write(var%lupri,'(/a,i4)') 'Print xkappa matrix, IZ =',iz
      do i=1,nloc
        write(var%lupri,'(i4,10es9.1)') i,(xkappa(i,j,iz),j=1,nloc)
      enddo
    enddo
  endif

  ! get new MO coefficients(scoef)
  if(restore)then
    scoef(:,:,:) = storec(:,:,:)  ! restore MO coefficients
  elseif(rejected)then
    scoef(:,:,:) = storec(:,:,:)  ! restore MO coefficients; the step was rejected
  elseif(.not.rejected)then
    storec(:,:,:) = scoef(:,:,:)  ! store MO coefficients
  endif
  call rotcmo_mod(var%nfbas(ifrp,0),nloc,4,xkappa,nloc,var%ipqtoq(ifrp,0))

  ! deallocate local variable
  deallocate( xkappa )

  !------------------------------------------------------------------------------!
  end subroutine mv2point

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  real*8 function get_func()
  !==============================================================================!
  !                       Get value of your functional                           !
  !==============================================================================!
  implicit none

  ! Global variables
  ! real*8,  intent(inout) :: smat(*)
  ! real*8,  intent(  out) :: pmat(nloc,nloc,4,nucdep)

  ! Internal variables
  integer :: nloc, nucdep
  real*8  :: dinvmg

  !------------------------------------------------------------------------------!
  nloc   = size(pmat,1)
  nucdep = size(pmat,4)

  if(PM%prjloc)then ! localization with projection analysis
    call get_overlap(smat)
    call orblp1(var%nfbas(1,1),var%nfbas(1,2),nloc,nucdep,var%maxref,dinvmg)
  else               ! localization with Mulliken analysis
    call orblo1(var%nfbas(1,1),var%nfbas(1,2),nloc,nucdep,dinvmg)
  endif
  dinvmg4out = dinvmg   ! global variable
  get_func   = dinvmg

  !------------------------------------------------------------------------------!
  end function get_func

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine get_hessian(hmat,full,totdiag,totnond)
  !==============================================================================!
  !                                                                              !
  !  Construction of full hessian                                                !
  !                                                                              ! 
  !     HMAT         A = d^2 G / dK* dK     B = d^2 G / dK dK          K         !
  !  / A   B  \       / A_uu   A_ub \        / B_uu   B_ub \      / K(unbared) \ !
  !  \ B*  A* /       \ A_bu   A_bb /        \ B_bu   B_bb /      \ K(bared)   / !
  !                                                                              ! 
  !==============================================================================!
  implicit none

  ! Global variables
  ! real*8, intent(in) :: pmat(nloc,nloc,4,nucdep)

  ! External variables
  logical, intent(in ) :: full
  real*8,  intent(out) :: totdiag, totnond

  real*8,  intent(out), pointer :: hmat(:,:)

  ! Internal variables
  integer :: i, j, k, l, m, n, ij, kl, ier
  integer :: nucdep, nloc, ndimh, nloc12, lupri
  real*8  :: pij1, pij2, pkml, p2kml, p2lmk, buf, buf1, buf2

  integer, allocatable :: ijn(:,:), iju(:), ijl(:)

  real*8,  allocatable, dimension(:,:,:) :: p1, p2, p3, p4
  real*8,  allocatable, dimension(:,:)   :: pd

  !------------------------------------------------------------------------------!
  lupri = var%lupri

  ! Get some dimensions
  nucdep = size(pmat,4)
  nloc   = size(pmat,1)
  nloc12 = ( nloc*nloc - nloc ) / 2
  ndimh  = 4*nloc12

  !------------------------------------------------------------------------------!
  ! Calculate real representation of hessian and get the diagonal dominance info !
  !------------------------------------------------------------------------------!
  totdiag = d0
  totnond = d0

  if(.not.full)then  ! label 111

    if(associated(hmat)) deallocate( hmat )
    allocate( hmat(ndimh,1), stat=ier ); call errall(ier,'get_hessian','hmat')
    hmat(:,:) = d0

    ! unbared-unbared A^R + B^R contribution
    do m = 1,nucdep
      ij = 0
      do j=1,nloc
        do i=j+1,nloc
          ij = ij + 1
          pij1 = pmat(i,j,1,m)
          buf  = pmat(i,i,1,m) - pmat(j,j,1,m)
          hmat(ij,1) = hmat(ij,1) + d2 * ( d4*pij1*pij1 - buf*buf )
        enddo
      enddo
    enddo

    ! bared-bared A^R + B^R contribution
    do m = 1,nucdep
      ij = nloc12
      do j=1,nloc
        do i=j+1,nloc
          ij = ij + 1
          pij1 = pmat(i,j,3,m); pij2 = pmat(i,j,4,m)
          buf  = pmat(i,i,1,m) - pmat(j,j,1,m)
          hmat(ij,1) = hmat(ij,1) + d2 * ( d2*(pij1*pij1 + pij2*pij2) - buf*buf )
        enddo
      enddo
    enddo

    ! unbared-unbared A^R - B^R contribution
    do m = 1,nucdep
      ij = 2*nloc12
      do j=1,nloc
        do i=j+1,nloc
          ij = ij + 1
          pij2 = pmat(i,j,2,m)
          buf  = pmat(i,i,1,m) - pmat(j,j,1,m)
          hmat(ij,1) = hmat(ij,1) + d2 * ( d4*pij2*pij2 - buf*buf )
        enddo
      enddo
    enddo

    ! bared-bared A^R + B^R contribution
    do m = 1,nucdep
      ij = 3*nloc12
      do j=1,nloc
        do i=j+1,nloc
          ij = ij + 1
          pij1 = pmat(i,j,3,m); pij2 = pmat(i,j,4,m)
          buf  = pmat(i,i,1,m) - pmat(j,j,1,m)
          hmat(ij,1) = hmat(ij,1) + d2 * ( d2*(pij1*pij1 + pij2*pij2) - buf*buf )
        enddo
      enddo
    enddo

    do i=1,ndimh
      hmat(i,1) = - hmat(i,1) / dfloat(nloc) ! minus sign --> we are looking for the minimum
      totdiag = totdiag + dabs( hmat(i,1) )
    enddo

    if(PM%iprloc.ge.3)then
      write(lupri,'(/a,2f10.5)') 'Print diagonal (A^R + B^R)_uu'
      write(lupri,'(8d12.4)') (hmat(i,1),i=1,nloc12)
      write(lupri,'(/a,2f10.5)') 'Print diagonal (A^R + B^R)_bb'
      write(lupri,'(8d12.4)') (hmat(i,1),i=nloc12+1,2*nloc12)
      write(lupri,'(/a,2f10.5)') 'Print diagonal (A^R - B^R)_uu'
      write(lupri,'(8d12.4)') (hmat(i,1),i=2*nloc12+1,3*nloc12)
      write(lupri,'(/a,2f10.5)') 'Print diagonal (A^R - B^R)_bb'
      write(lupri,'(8d12.4)') (hmat(i,1),i=3*nloc12+1,4*nloc12)
    endif


  else ! label 111

    if(associated(hmat)) deallocate( hmat )
    allocate( hmat(ndimh,ndimh), stat=ier ); call errall(ier,'get_hessian','hmat')
    hmat(:,:) = d0

    allocate( ijn(nloc,nloc),        stat=ier ); call errall(ier,'minimization_drv','ijn ')
    allocate( iju(nloc),             stat=ier ); call errall(ier,'minimization_drv','iju ')
    allocate( ijl(nloc),             stat=ier ); call errall(ier,'minimization_drv','ijl ')
    allocate( pd (nloc,nucdep),      stat=ier ); call errall(ier,'minimization_drv','pd  ')
    allocate( p1 (nloc,nloc,nucdep), stat=ier ); call errall(ier,'minimization_drv','p1  ')
    allocate( p2 (nloc,nloc,nucdep), stat=ier ); call errall(ier,'minimization_drv','p2  ')
    allocate( p3 (nloc,nloc,nucdep), stat=ier ); call errall(ier,'minimization_drv','p3  ')
    allocate( p4 (nloc,nloc,nucdep), stat=ier ); call errall(ier,'minimization_drv','p4  ')

    ! prepare some arrays for effective construction of the hessian
    do m=1,nucdep
      do i=1,nloc
        pd(i,m) = pmat(i,i,1,m)
        do j=1,nloc
          p1(j,i,m) = pmat(j,i,1,m)
          p2(j,i,m) = pmat(j,i,2,m)
          p3(j,i,m) = pmat(j,i,3,m)
          p4(j,i,m) = pmat(j,i,4,m)
        enddo
      enddo
    enddo

    ijn(:,:) = 0
    kl = 0
    do l=1,nloc
      do k=l+1,nloc
        kl = kl + 1
        ijn(k,l) = kl
        ijn(l,k) = ijn(k,l)
      enddo
      iju(l) = l - 1
      ijl(l) = l + 1
    enddo

    !----------------------------------------------------------------------------!
    !                       Construction of the HESSIAN                          !
    !----------------------------------------------------------------------------!
    !                           kl = < 1 , nloc12 >                              !
    !----------------------------------------------------------------------------!
    do m = 1,nucdep
      kl = 0
      do l=1,nloc
        do k=l+1,nloc
          kl = kl + 1
          buf  = d4*p1(k,l,m)
          p2kml= d2*pd(k,m) - pd(l,m)
          p2lmk= d2*pd(l,m) - pd(k,m)
 
          ! (A^R + B^R)_uu = Real{ d^2 G / dK*_ij dK_kl } + Real{ d^2 G / dK*_ij dK*_kl }
          !                = Real{ (32) }                 + Real{ (36) }
          do j=1,iju(k)     ! i==k
           hmat(ijn(j,k),kl) = hmat(ijn(j,k),kl) + buf*p1(j,k,m) + (p2kml-pd(j,m))*p1(j,l,m)
          enddo

          do i=ijl(l),nloc  ! j==l
           hmat(ijn(i,l),kl) = hmat(ijn(i,l),kl) + buf*p1(i,l,m) + (p2lmk-pd(i,m))*p1(i,k,m)
          enddo

          do j=1,iju(l)     ! i==l
           hmat(ijn(j,l),kl) = hmat(ijn(j,l),kl) - buf*p1(j,l,m) - (p2lmk-pd(j,m))*p1(j,k,m)
          enddo

          do i=ijl(k),nloc  ! j==k
           hmat(ijn(i,k),kl) = hmat(ijn(i,k),kl) - buf*p1(i,k,m) - (p2kml-pd(i,m))*p1(i,l,m)
          enddo

        enddo
      enddo
    enddo

    !----------------------------------------------------------------------------!
    !                     kl = < nloc12 + 1 , 2*nloc12 >                         !
    !----------------------------------------------------------------------------!
    do m = 1,nucdep
      kl = nloc12
      do l=1,nloc
        do k=l+1,nloc
          kl = kl + 1
          buf1 = d2*p3(k,l,m)
          buf2 = d2*p4(k,l,m)
          buf  = d2*buf1
          pkml =    pd(k,m) - pd(l,m)
          p2kml= d2*pd(k,m) - pd(l,m)
          p2lmk= d2*pd(l,m) - pd(k,m)

          ! (A^R + B^R)_ub = Real{ d^2 G / dK*_ij dK_kl } + Real{ d^2 G / dK*_ij dK*_kl }
          !                = Real{ (33) }                 + Real{ (37) }
          do j=1,iju(k)     ! i==k
           hmat(ijn(j,k),kl) = hmat(ijn(j,k),kl) + buf*p1(j,k,m) + (p2kml-pd(j,m))*p3(j,l,m)
          enddo

          do i=ijl(l),nloc  ! j==l
           hmat(ijn(i,l),kl) = hmat(ijn(i,l),kl) + buf*p1(i,l,m) - (p2lmk-pd(i,m))*p3(i,k,m)
          enddo

          do j=1,iju(l)     ! i==l
           hmat(ijn(j,l),kl) = hmat(ijn(j,l),kl) - buf*p1(j,l,m) + (p2lmk-pd(j,m))*p3(j,k,m)
          enddo

          do i=ijl(k),nloc  ! j==k
           hmat(ijn(i,k),kl) = hmat(ijn(i,k),kl) - buf*p1(i,k,m) - (p2kml-pd(i,m))*p3(i,l,m)
          enddo

          ! (A^R + B^R)_bb = Real{ d^2 G / dK*_ij dK_kl } + Real{ d^2 G / dK*_ij dK*_kl }
          !                = Real{ (35) }                 + Real{ (39) }
          n = nloc12
          do j=1,iju(k)     ! i==k
           hmat(ijn(j,k)+n,kl) = hmat(ijn(j,k)+n,kl) &
                               - buf1*p3(j,k,m) - buf2*p4(j,k,m) + (p2kml-pd(j,m))*p1(j,l,m)
          enddo

          do i=ijl(l),nloc  ! j==l
           hmat(ijn(i,l)+n,kl) = hmat(ijn(i,l)+n,kl) &
                               + buf1*p3(i,l,m) + buf2*p4(i,l,m) + (p2lmk-pd(i,m))*p1(i,k,m)
          enddo

          do j=1,iju(l)     ! i==l
           hmat(ijn(j,l)+n,kl) = hmat(ijn(j,l)+n,kl) &
                               + buf1*p3(j,l,m) + buf2*p4(j,l,m) + (p2lmk-pd(j,m))*p1(j,k,m)
          enddo

          do i=ijl(k),nloc  ! j==k
           hmat(ijn(i,k)+n,kl) = hmat(ijn(i,k)+n,kl) &
                               - buf1*p3(i,k,m) - buf2*p4(i,k,m) + (p2kml-pd(i,m))*p1(i,l,m)
          enddo

        enddo
      enddo
    enddo

    !----------------------------------------------------------------------------!
    !                     kl = < 2*nloc12 + 1 , 3*nloc12 >                       !
    !----------------------------------------------------------------------------!
    do m = 1,nucdep
      kl = 2*nloc12
      do l=1,nloc
        do k=l+1,nloc
          kl = kl + 1
          buf  = d4*p2(k,l,m)
          pkml =    pd(k,m) - pd(l,m)
          p2kml= d2*pd(k,m) - pd(l,m)
          p2lmk= d2*pd(l,m) - pd(k,m)

          ! (B^I - A^I)_uu = Imag{ d^2 G / dK*_ij dK*_kl } - Imag{ d^2 G / dK*_ij dK_kl }
          !                = - Imag{ (36) }                - Imag{ (32) }
          do j=1,iju(k)     ! i==k
           hmat(ijn(j,k),kl) = hmat(ijn(j,k),kl) - buf*p1(j,k,m) - (p2kml-pd(j,m))*p2(j,l,m)
          enddo

          do i=ijl(l),nloc  ! j==l
           hmat(ijn(i,l),kl) = hmat(ijn(i,l),kl) - buf*p1(i,l,m) + (p2lmk-pd(i,m))*p2(i,k,m)
          enddo

          do j=1,iju(l)     ! i==l
           hmat(ijn(j,l),kl) = hmat(ijn(j,l),kl) + buf*p1(j,l,m) - (p2lmk-pd(j,m))*p2(j,k,m)
          enddo

          do i=ijl(k),nloc  ! j==k
           hmat(ijn(i,k),kl) = hmat(ijn(i,k),kl) + buf*p1(i,k,m) + (p2kml-pd(i,m))*p2(i,l,m)
          enddo

          ! (B^I - A^I)_bu = Imag{ d^2 G / dK*_ij dK*_kl } - Imag{ d^2 G / dK*_ij dK_kl }
          !                = - Imag{ (38) }                - Imag{ (34) }
          n = nloc12
          do j=1,iju(k)     ! i==k
           hmat(ijn(j,k)+n,kl) = hmat(ijn(j,k)+n,kl) + buf*p3(j,k,m) - (p2kml-pd(j,m))*p4(j,l,m)
          enddo

          do i=ijl(l),nloc  ! j==l
           hmat(ijn(i,l)+n,kl) = hmat(ijn(i,l)+n,kl) - buf*p3(i,l,m) - (p2lmk-pd(i,m))*p4(i,k,m)
          enddo

          do j=1,iju(l)     ! i==l
           hmat(ijn(j,l)+n,kl) = hmat(ijn(j,l)+n,kl) - buf*p3(j,l,m) - (p2lmk-pd(j,m))*p4(j,k,m)
          enddo

          do i=ijl(k),nloc  ! j==k
           hmat(ijn(i,k)+n,kl) = hmat(ijn(i,k)+n,kl) + buf*p3(i,k,m) - (p2kml-pd(i,m))*p4(i,l,m)
          enddo

          ! (A^R - B^R)_uu = Real{ d^2 G / dK*_ij dK_kl } - Real{ d^2 G / dK*_ij dK*_kl }
          !                = Real{ (32) }                 - Real{ (36) }
          n = 2*nloc12
          do j=1,iju(k)     ! i==k
           hmat(ijn(j,k)+n,kl) = hmat(ijn(j,k)+n,kl) - buf*p2(j,k,m) + (p2kml-pd(j,m))*p1(j,l,m)
          enddo

          do i=ijl(l),nloc  ! j==l
           hmat(ijn(i,l)+n,kl) = hmat(ijn(i,l)+n,kl) + buf*p2(i,l,m) + (p2lmk-pd(i,m))*p1(i,k,m)
          enddo

          do j=1,iju(l)     ! i==l
           hmat(ijn(j,l)+n,kl) = hmat(ijn(j,l)+n,kl) + buf*p2(j,l,m) + (p2lmk-pd(j,m))*p1(j,k,m)
          enddo

          do i=ijl(k),nloc  ! j==k
           hmat(ijn(i,k)+n,kl) = hmat(ijn(i,k)+n,kl) - buf*p2(i,k,m) + (p2kml-pd(i,m))*p1(i,l,m)
          enddo

        enddo
      enddo
    enddo

    !----------------------------------------------------------------------------!
    !                  kl = < 3*nloc12 + 1 , 4*nloc12 >                          !
    !----------------------------------------------------------------------------!
    do m = 1,nucdep
      kl = 3*nloc12
      do l=1,nloc
        do k=l+1,nloc
          kl = kl + 1
          buf1 = d2*p3(k,l,m)
          buf2 = d2*p4(k,l,m)
          buf  = d2*buf2
          pkml =    pd(k,m) - pd(l,m)
          p2kml= d2*pd(k,m) - pd(l,m)
          p2lmk= d2*pd(l,m) - pd(k,m)

          ! (B^I - A^I)_ub = Imag{ d^2 G / dK*_ij dK*_kl } - Imag{ d^2 G / dK*_ij dK_kl }
          !                = - Imag{ (37) }                - Imag{ (33) }
          do j=1,iju(k)     ! i==k
           hmat(ijn(j,k),kl) = hmat(ijn(j,k),kl) - buf*p1(j,k,m) - (p2kml-pd(j,m))*p4(j,l,m)
          enddo

          do i=ijl(l),nloc  ! j==l
           hmat(ijn(i,l),kl) = hmat(ijn(i,l),kl) - buf*p1(i,l,m) + (p2lmk-pd(i,m))*p4(i,k,m)
          enddo

          do j=1,iju(l)     ! i==l
           hmat(ijn(j,l),kl) = hmat(ijn(j,l),kl) + buf*p1(j,l,m) - (p2lmk-pd(j,m))*p4(j,k,m)
          enddo

          do i=ijl(k),nloc  ! j==k
           hmat(ijn(i,k),kl) = hmat(ijn(i,k),kl) + buf*p1(i,k,m) + (p2kml-pd(i,m))*p4(i,l,m)
          enddo

          ! (B^I - A^I)_bb = Imag{ d^2 G / dK*_ij dK*_kl } - Imag{ d^2 G / dK*_ij dK_kl }
          !                = - Imag{ (39) }                - Imag{ (35) }
          n = nloc12
          do j=1,iju(k)     ! i==k
           hmat(ijn(j,k)+n,kl) = hmat(ijn(j,k)+n,kl) &
                               - buf1*p4(j,k,m) + buf2*p3(j,k,m) + (p2kml-pd(j,m))*p2(j,l,m)
          enddo

          do i=ijl(l),nloc  ! j==l
           hmat(ijn(i,l)+n,kl) = hmat(ijn(i,l)+n,kl) &
                               + buf1*p4(i,l,m) - buf2*p3(i,l,m) + (p2lmk-pd(i,m))*p2(i,k,m)
          enddo

          do j=1,iju(l)     ! i==l
           hmat(ijn(j,l)+n,kl) = hmat(ijn(j,l)+n,kl) &
                               + buf1*p4(j,l,m) - buf2*p3(j,l,m) + (p2lmk-pd(j,m))*p2(j,k,m)
          enddo

          do i=ijl(k),nloc  ! j==k
           hmat(ijn(i,k)+n,kl) = hmat(ijn(i,k)+n,kl) &
                               - buf1*p4(i,k,m) + buf2*p3(i,k,m) + (p2kml-pd(i,m))*p2(i,l,m)
          enddo

          ! (A^R - B^R)_ub = Real{ d^2 G / dK*_ij dK_kl } - Real{ d^2 G / dK*_ij dK*_kl }
          !                = Real{ (33) }                 - Real{ (37) }
          n = 2*nloc12
          do j=1,iju(k)     ! i==k
           hmat(ijn(j,k)+n,kl) = hmat(ijn(j,k)+n,kl) - buf*p2(j,k,m) + (p2kml-pd(j,m))*p3(j,l,m)
          enddo

          do i=ijl(l),nloc  ! j==l
           hmat(ijn(i,l)+n,kl) = hmat(ijn(i,l)+n,kl) + buf*p2(i,l,m) + (p2lmk-pd(i,m))*p3(i,k,m)
          enddo

          do j=1,iju(l)     ! i==l
           hmat(ijn(j,l)+n,kl) = hmat(ijn(j,l)+n,kl) + buf*p2(j,l,m) + (p2lmk-pd(j,m))*p3(j,k,m)
          enddo

          do i=ijl(k),nloc  ! j==k
           hmat(ijn(i,k)+n,kl) = hmat(ijn(i,k)+n,kl) - buf*p2(i,k,m) + (p2kml-pd(i,m))*p3(i,l,m)
          enddo

          ! (A^R - B^R)_bb = Real{ d^2 G / dK*_ij dK_kl } - Real{ d^2 G / dK*_ij dK*_kl }
          !                = Real{ (35) }                 - Real{ (39) }
          n = 3*nloc12
          do j=1,iju(k)     ! i==k
           hmat(ijn(j,k)+n,kl) = hmat(ijn(j,k)+n,kl) &
                               - buf1*p3(j,k,m) - buf2*p4(j,k,m) + (p2kml-pd(j,m))*p1(j,l,m)
          enddo

          do i=ijl(l),nloc  ! j==l
           hmat(ijn(i,l)+n,kl) = hmat(ijn(i,l)+n,kl) &
                               + buf1*p3(i,l,m) + buf2*p4(i,l,m) + (p2lmk-pd(i,m))*p1(i,k,m)
          enddo

          do j=1,iju(l)     ! i==l
           hmat(ijn(j,l)+n,kl) = hmat(ijn(j,l)+n,kl) &
                               + buf1*p3(j,l,m) + buf2*p4(j,l,m) + (p2lmk-pd(j,m))*p1(j,k,m)
          enddo

          do i=ijl(k),nloc  ! j==k
           hmat(ijn(i,k)+n,kl) = hmat(ijn(i,k)+n,kl) &
                               - buf1*p3(i,k,m) - buf2*p4(i,k,m) + (p2kml-pd(i,m))*p1(i,l,m)
          enddo

        enddo
      enddo
    enddo

    !----------------------------------------------------------------------------!
    !            Calculate diagonal dominance and deallocate arrays              !
    !----------------------------------------------------------------------------!
    do j=1,ndimh
      do i=1,j
        hmat(i,j) = - hmat(i,j) / dfloat(nloc)  ! minus sign --> we are looking for the minimum
        if(i==j)then
          totdiag = totdiag + dabs( hmat(i,j) )
        else
          totnond = totnond + dabs( hmat(i,j) )
        endif
      enddo
    enddo
    totnond = d2*totnond

    !----------------------------------------------------------------------------!
    !                              PRINT section                                 !
    !----------------------------------------------------------------------------!
    if(PM%iprloc.ge.3)then
      write(lupri,'(/a,2f10.5)') 'Print (A^R + B^R)_uu'
      do i=1,nloc12
        write(lupri,'(i4,15f10.5)') i,(hmat(i,j),j=1,nloc12)
      enddo

      write(lupri,'(/a,2f10.5)') 'Print (B^I - A^I)_uu'
      do i=1,nloc12
        write(lupri,'(i4,15f10.5)') i,(hmat(i,j),j=2*nloc12+1,3*nloc12)
      enddo

      write(lupri,'(/a,2f10.5)') 'Print (A^R - B^R)_uu'
      do i=2*nloc12+1,3*nloc12
        write(lupri,'(i4,15f10.5)') i,(hmat(i,j),j=2*nloc12+1,3*nloc12)
      enddo

      write(lupri,'(/a,2f10.5)') 'Print (A^R + B^R)_ub'
      do i=1,nloc12
        write(lupri,'(i4,15f10.5)') i,(hmat(i,j),j=nloc12+1,2*nloc12)
      enddo

      write(lupri,'(/a,2f10.5)') 'Print (B^I - A^I)_ub'
      do i=1,nloc12
        write(lupri,'(i4,15f10.5)') i,(hmat(i,j),j=3*nloc12+1,4*nloc12)
      enddo

      write(lupri,'(/a,2f10.5)') 'Print (A^R - B^R)_ub'
      do i=2*nloc12+1,3*nloc12
        write(lupri,'(i4,15f10.5)') i,(hmat(i,j),j=3*nloc12+1,4*nloc12)
      enddo

      write(lupri,'(/a,2f10.5)') 'Print (A^R + B^R)_bb'
      do i=nloc12+1,2*nloc12
        write(lupri,'(i4,15f10.5)') i,(hmat(i,j),j=nloc12+1,2*nloc12)
      enddo

      write(lupri,'(/a,2f10.5)') 'Print (B^I - A^I)_bb'
      do i=nloc12+1,2*nloc12
        write(lupri,'(i4,15f10.5)') i,(hmat(i,j),j=3*nloc12+1,4*nloc12)
      enddo

      write(lupri,'(/a,2f10.5)') 'Print (A^R - B^R)_bb'
      do i=3*nloc12+1,4*nloc12
        write(lupri,'(i4,15f10.5)') i,(hmat(i,j),j=3*nloc12+1,4*nloc12)
      enddo

      write(lupri,'(/a,2f10.5)') 'Print (B^I - A^I)_bu'
      do i=nloc12+1,2*nloc12
        write(lupri,'(i4,15f10.5)') i,(hmat(i,j),j=2*nloc12+1,3*nloc12)
      enddo
    endif

    ! deallocate arrays
    deallocate( ijn )
    deallocate( iju )
    deallocate( ijl )
    deallocate( pd  )
    deallocate( p1  )
    deallocate( p2  )
    deallocate( p3  )
    deallocate( p4  )

  endif ! label 111

  ! Calculate mean value
  totdiag = totdiag / dfloat(ndimh)
  totnond = totnond / dfloat(ndimh*ndimh-ndimh)

  end subroutine get_hessian

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine get_gradient(gradient,dnrmg)
  !==============================================================================!
  !                     Construction of the gradient                             !
  !==============================================================================!
  implicit none

  ! Global variables
  ! real*8, intent(in) :: pmat(nloc,nloc,4,nucdep)

  ! External variables
  real*8, intent(out) :: dnrmg

  real*8, intent(out), dimension(:) :: gradient

  ! Internal variables
  integer :: i, j, k, ij1, ij2, ij3, ij4
  integer :: nucdep, nloc, ndimg, nloc12, lupri

  real*8  :: dtmp, dnrm2

  !------------------------------------------------------------------------------!
  lupri = var%lupri
  gradient(:) = d0

  ! Get some dimensions
  nucdep = size(pmat,4)
  nloc   = size(pmat,1)
  ndimg  = size(gradient,1)
  nloc12 = ndimg / 4

  !------------------------------------------------------------------------------!
  !               Calculate real representation of the gradient                  !
  !------------------------------------------------------------------------------!
  do k = 1,nucdep
    ij1 = 0
    do j=1,nloc
      do i=j+1,nloc
        ij1 = ij1 + 1
        ij2 = ij1 + nloc12
        ij3 = ij2 + nloc12
        ij4 = ij3 + nloc12
        dtmp = d2*(pmat(j,j,1,k) - pmat(i,i,1,k))
        gradient(ij1) = gradient(ij1) + dtmp * pmat(i,j,1,k)
        gradient(ij2) = gradient(ij2) + dtmp * pmat(i,j,3,k)
        gradient(ij3) = gradient(ij3) - dtmp * pmat(i,j,2,k)
        gradient(ij4) = gradient(ij4) - dtmp * pmat(i,j,4,k)
      enddo
    enddo
  enddo

  do i=1,ndimg
    gradient(i) = - gradient(i) / dfloat(nloc)  ! minus sign --> we are looking for the minimum
  enddo

  if(PM%iprloc.ge.2)then
    write(lupri,'(/a,2f10.5)') 'Print gradient^R_u'
    write(lupri,'(8d12.4)') (gradient(i),i=1,nloc12)
    write(lupri,'(/a,2f10.5)') 'Print gradient^R_b}'
    write(lupri,'(8d12.4)') (gradient(i),i=1+nloc12,2*nloc12)
    write(lupri,'(/a,2f10.5)') 'Print gradient^I_u'
    write(lupri,'(8d12.4)') (gradient(i),i=1+2*nloc12,3*nloc12)
    write(lupri,'(/a,2f10.5)') 'Print gradient^I_b'
    write(lupri,'(8d12.4)') (gradient(i),i=1+3*nloc12,4*nloc12)
  endif

  ! Norm of the gradient
  dnrmg = dnrm2(ndimg,gradient,1)

  end subroutine get_gradient

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine get_date(daytid)
  implicit none
  character, intent(out) :: daytid*24

  call gtinfo(daytid)

  end subroutine get_date
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  SUBROUTINE ORBLP1(NBL,NBS,NLOC,NUCDEP,MAXREF,DINVM)
!***********************************************************************
!
!     Based on T.Saue March 26 1999 : SUBROUTINE ANAPRJ
!     Modify by S. Dubillard - Apr 21 2004
!
!***********************************************************************
      USE SELECT_COEFF_OWN_BAS
      IMPLICIT NONE

      ! Global variables
      ! REAL*8,  INTENT(  OUT) :: PMAT(NLOC,NLOC,4,NUCDEP)

      ! External variables
      INTEGER, INTENT(IN ) :: NBL, NBS, NLOC, NUCDEP, MAXREF
      REAL*8,  INTENT(OUT) :: DINVM

      ! Internal variables
      CHARACTER*6  :: REFFIL(MAXREF)
      CHARACTER*72 :: VECREF(2,MAXREF)
      REAL*8, PARAMETER :: DM1=-1.0D0

      LOGICAL   :: TOBE

      INTEGER :: I, J, K, IZ, IER, IFRP, ICFS, NB, NQ, NCDIM, NBDIM, NMOL
      INTEGER :: IVEC, JOFF, NREF, N2REF, LUPRI, NPSH, NESH, NORB, NREFS
      INTEGER :: NSTR(2,0:2,0:MAXREF), NVECS(2,0:2), IOFF(2), IDUMARR(2), IPQTOQ(4,0:7)

      REAL*8    :: DDUMARR(2), DINV, QLOC2, DLOC

      INTEGER, ALLOCATABLE, DIMENSION(:) :: KVECTMP, JREF, IPIV
      REAL*8,  ALLOCATABLE, DIMENSION(:) :: RCOEF, CMO, BV, GVEC, AMAT, WMAT

      TYPE :: PTR2ARRAY
        INTEGER, DIMENSION(:), POINTER :: P
      END TYPE PTR2ARRAY

      TYPE(PTR2ARRAY), ALLOCATABLE :: KVEC(:,:)

  !------------------------------------------------------------------------------!
      IFRP   = 1
      LUPRI= VAR%LUPRI
      NPSH = VAR%NPSH(IFRP)
      NESH = VAR%NESH(IFRP)
      NORB = VAR%NORB(IFRP)
      NREFS= VAR%NREFS
      IPQTOQ(:,:) = VAR%IPQTOQ(:,:)
      VECREF(:,:) = VAR_VECREF(:,:)
      REFFIL(:)   = VAR_REFFIL(:)

      ALLOCATE( KVEC(2,NREFS), STAT=IER )
      IF(IER /= 0)THEN
        WRITE(LUPRI,*) "Error in allocation of KVEC array of pointers in OBLP1"
        WRITE(LUPRI,*) "STAT=",IER
        WRITE(0,*) "Error in allocation ==> see output file "
        STOP
      ENDIF
      DO I=1,2
      DO J=1,NREFS
        NULLIFY(KVEC(I,J)%P)
      ENDDO
      ENDDO

      NB     = NBL + NBS
!     DO IFRP = 1,NFSYM
        NVECS(IFRP,1) = 0
        ALLOCATE( KVECTMP(NORB), STAT=IER ); CALL ERRALL(IER,'ORBLP1','KVECTMP')
        DO I = 1,NREFS
          NSTR(IFRP,0,I) = 1
          CALL NUMLST_MOD(VECREF(IFRP,I),KVECTMP, &
               NORB,-NPSH,NESH,               &
               IFRP,NSTR(IFRP,0,I))
          CALL ORBCNT(KVECTMP,NSTR(IFRP,0,I), &
               NPSH,NESH,                     &
               NSTR(IFRP,2,I),NSTR(IFRP,1,I))
          NSTR(IFRP,0,I) = NSTR(IFRP,1,I) + NSTR(IFRP,2,I)
          NVECS(IFRP,1) = NVECS(IFRP,1) + NSTR(IFRP,0,I)
          ALLOCATE(KVEC(IFRP,I)%P(NSTR(IFRP,0,I)),STAT=IER); CALL ERRALL(IER,'ORBLP1','KVEC%P')
          KVEC(IFRP,I)%P(1:NSTR(IFRP,0,I)) = KVECTMP(1:NSTR(IFRP,0,I))
        ENDDO
        DEALLOCATE( KVECTMP )
!     ENDDO
      NCDIM = NB*NVECS(IFRP,1)*4
      ALLOCATE( RCOEF(NCDIM), STAT=IER ); CALL ERRALL(IER,'ORBLP1','RCOEF')

!     COEFFICIENT SELECTION
!     =====================

      ICFS = 6
      IF(VAR%IPRPRJ.GE.5) ICFS = 7
      ALLOCATE( CMO(VAR%N2BBASXQ), STAT=IER ); CALL ERRALL(IER,'ORBLP1','CMO')

      IOFF(1) = 0
      VAR%IPRPRJ = 4
      DO I = 1,NREFS
        INQUIRE(FILE=REFFIL(I),EXIST=TOBE)
        IF (.NOT.TOBE) THEN
          WRITE(LUPRI,'(A,A6,A)') 'Reference coefficient file ',REFFIL(I),' not present'
          CALL QUIT('PRJANA: No reference coefficients !')
        ENDIF
        IF(VAR%OWNBAS)THEN
          CALL SELOWN_DRV(2,CMO,DDUMARR,IDUMARR,               &
               I,1,VAR%LUCOEF,REFFIL(I),RCOEF,DDUMARR,IDUMARR, &
               NVECS(1,1),IOFF,KVEC(IFRP,I)%P,                 &
               NSTR(1,0,I))
          IF(PM%IPRLOC.GE.5)THEN
            WRITE(LUPRI,'(A,A6)') '* Selected coefficients for reference ',REFFIL(I)
            NQ = NB*IOFF(IFRP)+1
            CALL PRQMAT(RCOEF(NQ),NB,NSTR(IFRP,0,I),NB,NSTR(IFRP,0,I),4,IPQTOQ(1,0),LUPRI)
          ENDIF
!         DO IFRP = 1, NFSYM
            IOFF(IFRP) = IOFF(IFRP) + NSTR(IFRP,0,I)
!         ENDDO
        ELSE
          WRITE(LUPRI,*) 'Current version work only with .OWNBAS option'
          WRITE(0,*) 'Current version work only with .OWNBAS option'
          STOP
        ENDIF
      ENDDO

!     Throw away the full coefficients and eigenvalues
!     ================================================

      DEALLOCATE( CMO )

!     Do projection analysis
!     ======================

!     NBDIM = 0
!     DO IFRP = 1,NFSYM
        NBDIM = NVECS(IFRP,1)*NLOC*4
!     ENDDO
      ALLOCATE( BV(NBDIM), STAT=IER ); CALL ERRALL(IER,'ORBLP1','BV')
!     DO IFRP = 1,NFSYM
!     CALL HEADER('Fermion ircop '//FREP(IFRP),-1)
      NMOL   = NLOC
      NREF   = NVECS(IFRP,1)
      N2REF  = NREF*NREF
      ALLOCATE( GVEC((NREF+1)*2), STAT=IER ); CALL ERRALL(IER,'ORBLP1','GVEC')
      ALLOCATE( AMAT(N2REF*4),    STAT=IER ); CALL ERRALL(IER,'ORBLP1','AMAT')
      ALLOCATE( JREF(NREF*3),     STAT=IER ); CALL ERRALL(IER,'ORBLP1','JREF')
      ALLOCATE( IPIV(NREF),       STAT=IER ); CALL ERRALL(IER,'ORBLP1','IPIV')
      IVEC = 0
      DO I = 1,NREFS
        CALL PRJINF(JREF,IFRP,I,IVEC,NSTR,KVEC(IFRP,I)%P)
        IVEC = IVEC + NSTR(IFRP,0,I)
      ENDDO
      CALL ORBLP2(RCOEF,AMAT,NB,NLOC,NVECS(IFRP,1),BV,NSTR,IPIV,MAXREF)

! Memory allocation for the matrix W = Ac

      ALLOCATE( WMAT(NREF*NLOC*4), STAT=IER ); CALL ERRALL(IER,'ORBLP1','WMAT')

! Construction of W = Ac

      IF(PM%IPRLOC.GE.5)THEN
        WRITE(LUPRI,'(/A)') '* Matrice BV'
        CALL PRQMAT(BV,NREF,NLOC,NREF,NLOC,4,IPQTOQ(1,0),LUPRI)
        WRITE(LUPRI,'(/A)') '* Matrice AMAT'
        CALL PRQMAT(AMAT,NREF,NREF,NREF,NLOC,4,IPQTOQ(1,0),LUPRI)
      ENDIF
      CALL QGEMM(NREF,NLOC,NREF,D1,              &
           'N','N',IPQTOQ(1,0),AMAT,NREF,NREF,4, &
           'N','N',IPQTOQ(1,0),BV,NREF,NLOC,4,   &
           D0,IPQTOQ(1,0),WMAT,NREF,NLOC,4)

      IF(PM%IPRLOC.GE.5)THEN
        WRITE(LUPRI,'(/A)') '* Matrice W'
        CALL PRQMAT(WMAT,NREF,NLOC,NREF,NLOC,4,IPQTOQ(1,0),LUPRI)
      ENDIF

! Construct P_A for each center A

      JOFF = 0
      DO I = 1,NUCDEP
        CALL QGEMM(NLOC,NLOC,NSTR(IFRP,0,I),D1,                  &
                   'H','N',IPQTOQ(1,0),BV(JOFF+1),NREF,NLOC,4,   &
                   'N','N',IPQTOQ(1,0),WMAT(JOFF+1),NREF,NLOC,4, &
                   D0,IPQTOQ(1,0),PMAT(1,1,1,I),NLOC,NLOC,4)
        JOFF = JOFF + NSTR(IFRP,0,I)
      ENDDO

! P_A matrices symmetrization

      DO I = 1,NUCDEP
        CALL FULMAT('S',NLOC,NLOC,PMAT(1,1,1,I))
        DO IZ = 2,4
          CALL FULMAT('A',NLOC,NLOC,PMAT(1,1,IZ,I))
        ENDDO
      ENDDO

! Multiply pmat with 1/2 just to avoid confusion in the code below
      DO IZ=1,4
        DO I=1,NLOC
          DO J=1,NLOC
            DO K=1,NUCDEP
              PMAT(I,J,IZ,K) = PMAT(I,J,IZ,K) / D2
            ENDDO
          ENDDO
        ENDDO
      ENDDO

      IF(PM%IPRLOC.GE.2)THEN
        DO IZ=1,4
          WRITE(LUPRI,*)
          WRITE(LUPRI,'(A,I3)') 'Print Q^A matrices, IZ = ',IZ
          DO I=1, NUCDEP
            WRITE(LUPRI,'(A,I4)') ' Nucleus #',I
            DO J=1,NLOC
              WRITE(LUPRI,'(I4,10ES9.1)') J,(PMAT(J,K,IZ,I),K=1,NLOC)
            ENDDO
            WRITE(LUPRI,*)
          ENDDO
        ENDDO
      ENDIF

      IF(PM%IPRLOC.GE.5)THEN
        DO I = 1,NUCDEP
          WRITE(LUPRI,'(/A,I3)') '* ORBLP1: Matrice P  Atom ', I
          WRITE(LUPRI,'(A)') '  ============================'
          CALL PRQMAT(PMAT(1,1,1,I),NLOC,NLOC,NLOC,NLOC,4,IPQTOQ(1,0),LUPRI)
        ENDDO
      ENDIF

! Calculation and print of the localisation criterion

    WRITE(LUPRI,'(/20X,A)') '** Pipek-Mezey localization criteria **' 
    WRITE(LUPRI,'(20X,A)')  '  ===================================  '
    DINV = D0
    DO I = 1,NLOC
      QLOC2 = D0
      DO J = 1,NUCDEP
         QLOC2 = QLOC2 + PMAT(I,I,1,J)*PMAT(I,I,1,J)
      ENDDO
      DLOC = D1/QLOC2
      DINV = DINV + QLOC2
      WRITE(LUPRI,'(11X,A,I4,A,F7.4)')' Orbital',I,' Pipek-Mezey localization measure d = ',DLOC
    ENDDO
    DINVM = DINV/DFLOAT(NLOC)
    WRITE(LUPRI,'(/11X,A,F21.15)')' Inverse mean delocalization = ',DINVM
    DINVM = - DINVM  ! we are looking for the minimum

  !------------------------------------------------------------------------------!
  !                             Deallocate arrays                                !
  !------------------------------------------------------------------------------!
    DO I=1,NREFS
      DEALLOCATE( KVEC(1,I)%P )
    ENDDO
    DEALLOCATE( KVEC  )
    DEALLOCATE( RCOEF )
    DEALLOCATE( BV    )
    DEALLOCATE( GVEC  )
    DEALLOCATE( AMAT  )
    DEALLOCATE( JREF  )
    DEALLOCATE( IPIV  )
    DEALLOCATE( WMAT  )

  END SUBROUTINE ORBLP1

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  SUBROUTINE ORBLP2(RCOEF,AMAT,NB,NMOL,NREF,BVEC,NSTR,IPIVOT,MAXREF)
!***********************************************************************
!
!     Projection analysis for fermion ircop IFRP
!
!     Written by T.Saue March 26 1999
!
!***********************************************************************
      IMPLICIT NONE

      ! Global variables
      ! REAL*8,  INTENT(IN   ) :: SCOEF(:,:,:)
      ! REAL*8,  INTENT(INOUT) :: SMAT(*)

      ! External variables
      INTEGER, INTENT(IN ) :: NB, NMOL, NREF, MAXREF
      INTEGER, INTENT(IN ) :: NSTR(2,0:2,0:MAXREF)
      INTEGER, INTENT(OUT) :: IPIVOT(*)

      REAL*8,  INTENT(IN   ) :: RCOEF(*)
      REAL*8,  INTENT(  OUT) :: AMAT(NREF,NREF,4), BVEC(NREF,NMOL,4)

      ! Internal variables
      CHARACTER*6  :: REFFIL(MAXREF)
      CHARACTER*72 :: VECREF(2,MAXREF)
      INTEGER      :: I, J, IFRP, JOFF, KOFF, JOB, NEFF, NREFS, LUPRI
      INTEGER      :: IPQTOQ(4,0:7)

      REAL*8 :: PTOL
      REAL*8, PARAMETER :: DC=1.0D-2

      LOGICAL DOPOP

  !------------------------------------------------------------------------------!
      IFRP = 1
      NREFS = VAR%NREFS

      DOPOP = NREFS.GT.1
      PTOL  = VAR%PROTHR*DC
      LUPRI = VAR%LUPRI
      IPQTOQ(:,:) = VAR%IPQTOQ(:,:)
      VECREF(:,:) = VAR_VECREF(:,:)
      REFFIL(:)   = VAR_REFFIL(:)

!     Make projection vector

      CALL QTRANS90_MOD('AOMO','S',D0,NB,NB,NREF,NMOL, &
                        SMAT,NB,NB,1,IPQTOQ(1,0),      &
                        BVEC,NREF,NMOL,4,IPQTOQ(1,0),  &
                        RCOEF,NB,NREF,4,IPQTOQ(1,0),   &
                        SCOEF,NB,NMOL,4,IPQTOQ(1,0),VAR%IPRPRJ)
      IF(NREFS.GT.1) THEN

!       Make overlap matrix

        CALL QTRANS90_MOD('AOMO','S',D0,NB,NB,NREF,NREF, &
                          SMAT,NB,NB,1,IPQTOQ(1,0),      &
                          AMAT,NREF,NREF,4,IPQTOQ(1,0),  &
                          RCOEF,NB,NREF,4,IPQTOQ(1,0),   &
                          RCOEF,NB,NREF,4,IPQTOQ(1,0),VAR%IPRPRJ)

!       Print overlaps

        IF(PM%IPRLOC.GE.5)THEN
           JOFF = 1
           DO J = 1,NREFS
              KOFF = 1
              DO I = 1,(J-1)
                 WRITE(LUPRI,'(A,2(3X,A6))') '*Overlaps :',REFFIL(I),REFFIL(J)
                 CALL PRQMAT(AMAT(KOFF,JOFF,1),                  &
                      NSTR(IFRP,0,I),NSTR(IFRP,0,J),NREF,NREF,4, &
                      IPQTOQ(1,0),LUPRI)
                 KOFF = KOFF + NSTR(IFRP,0,I)
              ENDDO
              JOFF = JOFF + NSTR(IFRP,0,J)
           ENDDO
        ENDIF

!       Solve linear system by Cholesky decomposition
!       with full pivoting
!        --> notice that overlap matrix SMAT is destroyed !

        JOB = 1
        CALL QCHOLD(AMAT,NREF,4,NREF,NREF,SMAT,PTOL,NEFF,JOB,IPIVOT)
        CALL QCHOLS(AMAT,NREF,NEFF,NMOL,4,NREF,NREF,SMAT,BVEC, &
             NREF,NMOL,BVEC,NREF,NMOL,JOB,IPIVOT,SMAT(1+NREF))
        IF(NEFF.LT.NREF) THEN
          WRITE(LUPRI,'(A,A)') '* WARNING: ', &
           ' linear dependencies detected in Cholesky decomposition'
          WRITE(LUPRI,'(11X,A,I5,A,I5)') ' Reduced dimensionality: ',NREF,' --> ',NEFF
        ENDIF

!       Restore overlap matrix

        CALL QHMRST(AMAT,NREF,VAR%NZ,NREF,NREF)
      ENDIF

      END SUBROUTINE ORBLP2

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine orblo1(nbl,nbs,nloc,nucdep,dinvm)
  !==============================================================================!
  !                                                                              !
  !    Construction of Q^A_i matrices from Mulliken analysis (eq. 3 in [1])      !
  !                         Q^A_i == PMAT(i,i,4,A) =                             !
  !     = 1/2 ( SUM_{nu=A} SUM_{mu} C_{nu,i}^* S_{nu,mu} C_{mu,i} + c.c. )       !
  !                                                                              !
  !              Calculation of the inverse mean delocalization                  !
  !                  DINVM = 1/M  SUM_i^M SUM_A ( Q^A_i )^2                      !
  !==============================================================================!
  implicit none

  ! Global variables
  ! real*8,  intent(in ) :: smat(*), scoef(:,:,:)
  ! real*8,  intent(out) :: pmat(nloc,nloc,4,nucdep)
  ! integer, intent(in ) :: nfnuc(*)

  ! External variables
  integer, intent(in ) :: nbl, nbs, nloc, nucdep
  real*8,  intent(out) :: dinvm

  ! Internal variables
  integer :: i, j, k, iz, nb, ioffl, ioffs, lupri, ier
  integer :: ipqtoq(4,0:7)
  real*8  :: dinv, qloc2, dloc

  real*8, allocatable, dimension(:) :: wlmat, wsmat

  !------------------------------------------------------------------------------!
  nb = nbl + nbs

  lupri       = var%lupri
  ipqtoq(:,:) = var%ipqtoq(:,:)

  ! Memory allocation for the matrix W = Sc
  allocate( wlmat(nbl*nloc*4), stat=ier ); call errall(ier,'orblo1','wlmat')
  allocate( wsmat(nbs*nloc*4), stat=ier ); call errall(ier,'orblo1','wsmat')

  ! Construction of W = Sc (large and small components)
  call qgemm(nbl,nloc,nbl,d1,                     &
             'N','N',ipqtoq(1,0),smat,nb,nb,1,    &
             'N','N',ipqtoq(1,0),scoef,nb,nloc,4, &
              d0,ipqtoq(1,0),wlmat,nbl,nloc,4)
  if(.not.var%h2comp)then
    call qgemm(nbs,nloc,nbs,d1,                                &
               'N','N',ipqtoq(1,0),smat(1+nb*nbl+nbl),nb,nb,1, &
               'N','N',ipqtoq(1,0),scoef(1+nbl,1,1),nb,nloc,4, &
                d0,ipqtoq(1,0),wsmat,nbs,nloc,4)
  endif

  if(pm%iprloc.ge.5)then
    write(lupri,'(/a)') '* WL matrix :'
    call prqmat(wlmat,nbl,nloc,nbl,nloc,4,ipqtoq(1,0),lupri)
    if(.not.var%h2comp)then
      write(lupri,'(/a)') '* WS matrix :'
      call prqmat(wsmat,nbs,nloc,nbs,nloc,4,ipqtoq(1,0),lupri)
    endif
  endif

  ! Construct P_A for each center A
  ioffl = 0
  ioffs = 0
  do i = 1,nucdep
    call qgemm(nloc,nloc,nfnuc(i),d1,                            &
               'H','N',ipqtoq(1,0),scoef(1+ioffl,1,1),nb,nloc,4, &
               'N','N',ipqtoq(1,0),wlmat(1+ioffl),nbl,nloc,4,    &
                d0,ipqtoq(1,0),pmat(1,1,1,i),nloc,nloc,4)
    if(.not.var%h2comp)then
      call qgemm(nloc,nloc,nfnuc(i+nucdep),d1,                         &
                 'H','N',ipqtoq(1,0),scoef(1+nbl+ioffs,1,1),nb,nloc,4, &
                 'N','N',ipqtoq(1,0),wsmat(1+ioffs),nbs,nloc,4,        &
                 d1,ipqtoq(1,0),pmat(1,1,1,i),nloc,nloc,4)
      ioffs = ioffs + nfnuc(nucdep+i)
    endif
    ioffl = ioffl + nfnuc(i)
  enddo

  ! P_A matrices symetrisation
  do i = 1,nucdep
    call fulmat('S',nloc,nloc,pmat(1,1,1,i))
    do iz = 2,4
      call fulmat('A',nloc,nloc,pmat(1,1,iz,i))
    enddo
  enddo

  ! Multiply pmat with 1/2 just to avoid confusion in the code below
  do iz=1,4
    do i=1,nloc
      do j=1,nloc
        do k=1,nucdep
          pmat(i,j,iz,k) = pmat(i,j,iz,k) / d2
        enddo
      enddo
    enddo
  enddo

  if(pm%iprloc.ge.2)then
    do iz=1,4
      write(lupri,'(/a,i3)') 'Print Q^A matrices, IZ = ',iz
      do i=1,nucdep
        write(lupri,'(a,i4)') ' Nucleus #',i
        do j=1,nloc
          write(lupri,'(i4,10es9.1)') j,(pmat(j,k,iz,i),k=1,nloc)
        enddo
        write(lupri,*)
      enddo
    enddo
  endif

  if(pm%iprloc.ge.5)then
     do i = 1,nucdep
        write(lupri,'(/a,i3)') '* ORBLO1: Matrice P  Atom ',i
        write(lupri,'(a)') '  ============================'
        call prqmat(pmat(1,1,1,i),nloc,nloc,nloc,nloc,4,ipqtoq(1,0),lupri)
     enddo
  endif

  ! Calculation and print of the localisation criterion
  write(lupri,'(/20x,a)') '** Pipek-Mezey localization criteria **' 
  write(lupri,'(20x,a)')  '  ===================================  '
  dinv = d0
  do i = 1,nloc
    qloc2 = d0
    do j = 1,nucdep
       qloc2 = qloc2 + pmat(i,i,1,j)*pmat(i,i,1,j)
    enddo
    dloc = d1/qloc2
    dinv = dinv + qloc2
    write(lupri,'(11x,a,i4,a,f7.4)')' Orbital',i,' Pipek-Mezey localization measure d = ',dloc
  enddo
  dinvm = dinv/dfloat(nloc)
  write(lupri,'(/11x,a,f21.15)')' Inverse mean delocalization = ',dinvm
  dinvm = - dinvm  ! we are looking for the minimum

  !------------------------------------------------------------------------------!
  !                             Deallocate arrays                                !
  !------------------------------------------------------------------------------!
  deallocate( wlmat )
  deallocate( wsmat )

  end subroutine orblo1

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  SUBROUTINE ORBLO4(INDNAO,NUCDEP)
!***********************************************************************
!
!     Written by S.Dubillard Jan 21 2004
!     Last revision:
!
!     Construction of the array NFNUC which contains the number of
!     large components basis functions for each center and then the 
!     number of small components basis functions for each center
!
!***********************************************************************
      IMPLICIT NONE

      ! Global variables
      ! INTEGER, INTENT(OUT) :: NFNUC(2*NUCDEP)

      ! External variables
      INTEGER, INTENT(IN ) :: NUCDEP
      INTEGER, INTENT(IN ) :: INDNAO(2,*)

      ! Internal variables
      INTEGER :: IFRP, I

  !------------------------------------------------------------------------------!
      IFRP = 1
      DO I = 1,NUCDEP-1
         NFNUC(I)=INDNAO(1,I+1)-INDNAO(1,I)
      ENDDO
      NFNUC(NUCDEP)=INDNAO(2,1)-INDNAO(1,NUCDEP)
      DO I = 1,NUCDEP-1
         NFNUC(I+NUCDEP)=INDNAO(2,I+1)-INDNAO(2,I)
      ENDDO
      NFNUC(2*NUCDEP)=VAR%NFBAS(IFRP,0)-INDNAO(2,NUCDEP)

  END SUBROUTINE ORBLO4

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine wcmo2file
  !==============================================================================!
  !                                                                              !
  !          Write the localised MO coefficients to the DFCOEF file              !
  !                                                                              !
  !==============================================================================!
   implicit none

   ! Global variables
   ! integer :: isel(:)      ! isel(nloc)
   ! real*8  :: scoef(:,:,:) ! scoef(nb,nloc,4)

   ! Internal variables
   integer :: ifrp, nz, no, nb, nloc, ier
   integer :: i, ii, iz, ioff,iopt,ioer
   real*8  :: toterg

   integer, allocatable :: ibeig(:)
   real*8,  allocatable :: cmo(:,:,:), eig(:)

  !------------------------------------------------------------------------------!
  ifrp = 1
  nz   = 4
  no   = var%norb(ifrp)
  nb   = size(scoef,1)    ! nfbas(ifrp,0)
  nloc = size(scoef,2)

  ! allocate arrays
  allocate( cmo  (nb,no,4),   stat=ier ); call errall(ier,'wcmo2file','cmo  ')
  allocate( eig  (var%norbt), stat=ier ); call errall(ier,'wcmo2file','eig  ')
  allocate( ibeig(var%norbt), stat=ier ); call errall(ier,'wcmo2file','ibeig')

  ! read EIG, IBEIG and TOTERG
  !  call reacmo_mod(var%lucoef,cmo,eig,ibeig,toterg,2)
  open(var%lucoef,file='DFCOEF',status='old',form='unformatted',access='sequential',iostat=ioer)
  if(ioer /= 0) call quit('REACMO_MOD: ERROR while opening DFCOEF or DFPCMO')
  iopt = 14
  call REACMO(var%lucoef,'DFCOEF',cmo,eig,ibeig,toterg,iopt)

  ! reorganize MO
  ioff = var%npsh(ifrp)
  if(var%h2comp .or. var%levyle) ioff = 0
  do iz = 1,nz
    do i = 1,nloc
      ii = ioff + isel(i)
      call dcopy(nb,scoef(1,i,iz),1,cmo(1,ii,iz),1)
    enddo
  enddo
  if (pm%iprloc.ge.5)then
    write(var%lupri,'(/a)') '*wcmo2file : New CMO matrix'
    call prqmat(cmo,nb,no,nb,no,4,var%ipqtoq(1,0),var%lupri)
  endif

  ! write MO to file and close him
  !  call wricmo_mod(var%lucoef,cmo,eig,ibeig,toterg)
  call WRICMO(var%lucoef,CMO,EIG,IBEIG,TOTERG,.FALSE.) ! SSYM=F in C1
  close(var%lucoef,status='keep')

  ! deallocate arrays
  deallocate( cmo   )
  deallocate( eig   )
  deallocate( ibeig )

  end subroutine wcmo2file

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine wricmo_mod(iunit,cmo,eig,ibeig,toterg)
  !==============================================================================!
  !                                                                              !
  !        Write DHF-coefficients and eigenvalues to (un)formatted file          !
  !             (based on wricmo routine from the main Dirac code)               !
  !                                                                              !
  !        idim(1,j) - number of positronic solutions                            !
  !        idim(2,j) - number of electronic solutions                            !
  !        idim(3,j) - number of AO-basis functions                              !
  !                                                                              !
  !==============================================================================!
  implicit none

  ! External variables
  integer, intent(in) :: iunit, ibeig(:)
  real*8,  intent(in) :: toterg, eig(:), cmo(:,:,:)

  ! Internal variables
  character :: text*74
  integer   :: i, j, ier, idim(3,2)

  !------------------------------------------------------------------------------!
  ! Open file
  open(iunit,file='DFCOEF',status='unknown',form='unformatted',access='sequential',iostat=ier)

  if(ier /= 0) call quit('WRICMO_MOD: ERROR while creating new DFCOEF file')
  rewind iunit

  ! write cmo, eig and ibeig to file
  do i = 1,var%nfsym
    idim(1,i) = var%npsh(i)
    idim(2,i) = var%nesh(i)
    idim(3,i) = var%nfbas(i,0)
  enddo
  text(1:50) = var%title
  call gtinfo(text(51:74))
  rewind iunit
  write(iunit) text,var%nfsym,((idim(i,j),i = 1,3),j=1,var%nfsym),toterg
  write(iunit) cmo
  write(iunit) eig
  write(iunit) ibeig

  ! Close file
  close(iunit,status='keep')

  end subroutine wricmo_mod

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine reacmo_mod(iunit,cmo,eig,ibeig,toterg,iopt)
  !==============================================================================!
  !                                                                              !
  !     Read DHF-coefficients and eigenvalues from unformatted file              !
  !                                                                              !
  !     Global variable:                                                         !
  !     keyform == 1  -->  unformatted input                                     !
  !     keyform == 2  -->  formatted input                                       !
  !                                                                              !
  !     IOPT = 1  -->  read only MO coefficients                                 !
  !     IOPT > 1  -->  read all information                                      !
  !                                                                              !
  !     IDIM(1,IFRP) - number of positronic solutions                            !
  !     IDIM(2,IFRP) - number of electronic solutions                            !
  !     IDIM(3,IFRP) - number of AO-basis functions                              !
  !                                                                              !
  !     Based on REACMO routine from the main Dirac code                 !
  !                                                                              !
  !==============================================================================!
  implicit none

  ! External variables
  integer,   intent( in) :: iunit, iopt
  integer,   intent(out) :: ibeig(:)
  real*8,    intent(out) :: cmo(var%ncmotq), eig(:)

  ! Internal variables
  character :: text*74
  logical   :: info
  integer   :: i, j, ifrp, ierr, ioer
  integer   :: nsym, nfsym, lupri
  integer   :: idim(3,2)
  real*8    :: toterg

  !------------------------------------------------------------------------------!
  ! set some parameters
  nfsym = var%nfsym
  lupri = var%lupri

  ! Open file
  ! =========
  if(keyform == 1)then
    open(iunit,file='DFCOEF',status='old',form='unformatted',access='sequential',iostat=ioer)
  elseif(keyform == 2)then
    open(iunit,file='DFPCMO',status='old',form='formatted',access='sequential',iostat=ioer)
  else
    write(var%lupri,'(a)') ' ERROR: wrong keyform keyword in reacmo_mod'
    stop
  endif

  if(ioer /= 0) call quit('REACMO_MOD: ERROR while opening DFCOEF or DFPCMO')
  rewind iunit

  ! Read title line
  ! ===============

  if(keyform == 1)then
    read(iunit,iostat=ioer) text,nsym,((idim(i,j),i = 1,3),j=1,nsym),toterg
  elseif(keyform == 2)then
    read(iunit,'(a74)',iostat=ioer) text
    read(iunit,*,iostat=ioer) nsym,((idim(i,j),i = 1,3),j=1,nsym)
    read(iunit,*,iostat=ioer) toterg
  endif
  if(ioer == -1) call quit('REACMO_MOD: END OF FILE reading TEXT')
  if(ioer /=  0) call quit('REACMO_MOD: ERROR reading TEXT')

  ! Compatibility check for NFSYM and NFBAS
  ! =======================================

  ierr = 0
  if(nsym .ne. nfsym) then
    ierr = ierr + 1
    info = .true.
    write(lupri,'(/a/a,2i5)') &
       'REACMO FATAL ERROR: Incompatible number of fermion ircops', &
       ' from file and in this run:',nsym,nfsym
  else
    do ifrp = 1,nfsym
      if(idim(3,ifrp).ne.var%nfbas(ifrp,0)) ierr = ierr + 1
    end do
    if (ierr.gt.0) then
      info = .true.
      write(lupri,'(/a/a,2(/a,i2,2i5))') &
      'REACMO FATAL ERROR. Incompatible number of basis functions', &
      ' from file and in this run:', (' Fermion ircop no.',ifrp,    &
      idim(3,ifrp),var%nfbas(ifrp,0),ifrp=1,nfsym)
    end if
  end if

  if (ierr.gt.0) call quit('REACMO: old DFCOEF file is not consistent with this calculation')

  ! Read coefficients
  ! =================
 
  if(keyform == 1) read(iunit,iostat=ioer) cmo
  if(keyform == 2) read(iunit,*,iostat=ioer) cmo

  if(ioer == -1) call quit('REACMO_MOD: END OF FILE reading coefficients')
  if(ioer /=  0) call quit('REACMO_MOD: ERROR reading coefficients')

  ! Read eigenvalues and single group irrep (spinfree and Levy-Leblond calcs)
  ! =========================================================================

  if(iopt > 1)then
    if(keyform == 1) read(iunit,iostat=ioer) eig
    if(keyform == 2) read(iunit,*,iostat=ioer) eig
    if(ioer == -1) call quit('REACMO_MOD: END OF FILE reading eigenvalues')
    if(ioer /=  0) call quit('REACMO_MOD: ERROR reading eigenvalues')

    if(keyform == 1) read(iunit,iostat=ioer) ibeig
    if(keyform == 2) read(iunit,*,iostat=ioer) ibeig
    if(ioer == -1) call quit('REACMO_MOD: END OF FILE reading irrep identification')
    if(ioer /=  0) call quit('REACMO_MOD: ERROR reading irrep identification')
  endif

  ! Close file
  ! ==========
  close(iunit,status='keep')

  end subroutine reacmo_mod

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  SUBROUTINE ROTCMO_MOD(NBASI,NORBI,NZ,XKAPPA,NORBT,IP)
!***********************************************************************
!
!     Rotate coefficients
!
!     SCOEF_{new} = SCOEF_{old} exp( -\kappa )
!     We use a trick a la Sostrup notes, chapter 3.1.5
!
!     Input:
!        SCOEF   - old coefficients
!        XKAPPA  - Kappa matrix (Kappa is overwritten with junk)
!        IPRINT  - print level
!
!     Output:
!        SCOEF   - new coefficients
!
!     Written by J. Thyssen - Nov 24 1998
!
!***********************************************************************
      IMPLICIT NONE

      ! Global variables
      ! REAL*8, INTENT(INOUT) :: SCOEF(:,:,:)  ! SCOEF(NBASI,NORBI,NZ)

      ! External variables
      INTEGER, INTENT(IN) :: NBASI, NORBI, NZ, NORBT, IP(4)
      REAL*8,  INTENT(IN) :: XKAPPA(NORBT,NORBT,NZ)

      ! Internal variables
      REAL*8, PARAMETER :: DTHRS = 1.0D-8, DSIXTH = D1/6.00D00, DP5 = 0.50D00
      REAL*8, PARAMETER :: D1P5 = 1.5D00, THREQL = 1.0D-12, DM1 = -D1

      CHARACTER SECTID*12,CPUTID*12,WALLTID*12

      INTEGER :: I, J, NORBXQ, IERR, IADDR, IOFF, I_MAX, J_MAX, LUPRI, IER
      REAL*8  :: CPU1, WALL1, CPU2, WALL2, CPU3, WALL3, CPU4, WALL4, EIG, SQEIG, D_MAX
      REAL*8  :: ABS, SQRT, COS, SIN

      REAL*8, ALLOCATABLE, DIMENSION(:) :: REIG, RSIN, RCOS, RV, RTMP

  !------------------------------------------------------------------------------!
      LUPRI  = VAR%LUPRI

      CALL GETTIM(CPU1,WALL1)
      IF (PM%IPRLOC.GE.5) CALL HEADER('Output from ROTCMO_MOD',-1)

      NORBXQ = NORBI*NORBI*NZ
      ALLOCATE( REIG(NORBI ), STAT=IER ); CALL ERRALL(IER,'ROTCMO_MOD','REIG')
      ALLOCATE( RSIN(NORBXQ), STAT=IER ); CALL ERRALL(IER,'ROTCMO_MOD','RSIN')
      ALLOCATE( RV  (NORBXQ), STAT=IER ); CALL ERRALL(IER,'ROTCMO_MOD','KV  ')
      ALLOCATE( RTMP(NORBXQ), STAT=IER ); CALL ERRALL(IER,'ROTCMO_MOD','RTMP')

      ALLOCATE( RCOS(NORBI*NBASI*NZ), STAT=IER ); CALL ERRALL(IER,'ROTCMO_MOD','RCOS')

!     Calculate -\kappa^2 = \kappa^{\dagger} \kappa
!     (use RCOS as temporary array)
!     ---------------------------------------------

      CALL QGEMM(NORBI,NORBI,NORBI,D1,             &
                 'H','N',IP,XKAPPA,NORBT,NORBT,NZ, &
                 'N','N',IP,XKAPPA,NORBT,NORBT,NZ, &
                 D0,IP,RCOS,NORBI,NORBI,NZ)

      IF (PM%IPRLOC.GE.10)THEN
         CALL HEADER('ROTCMO_MOD: -kappa^2',-1)
         CALL PRQMAT(RCOS,NORBI,NORBI,NORBI,NORBI,NZ,IP,LUPRI)
      END IF

!     Diagonalize -\kappa^2
!     ---------------------------------------------

      CALL QDIAG90_MOD(NZ,NORBI,RCOS,NORBI,NORBI, &
                       REIG,1,RV,NORBI,NORBI,IERR)
      IF (PM%IPRLOC.GE.10)THEN
         CALL HEADER('ROTCMO_MOD: V matrix (eigenvector)',-1)
         CALL PRQMAT(RV,NORBI,NORBI,NORBI,NORBI,NZ,IP,LUPRI)
         CALL HEADER('ROTCMO_MOD: eigenvalues',-1)
         CALL PRQMAT(REIG,1,NORBI,1,NORBI,1,IP,LUPRI)
      END IF

!     Make cos(d) and sin(d)/d ; d = sqrt(eigenvalues)
!     ------------------------------------------------

      CALL GETTIM(CPU2,WALL2)
      RCOS(1:NORBXQ) = D0
      RSIN(1:NORBXQ) = D0
      DO I = 1,NORBI
         IADDR = (I-1) + (I-1)*NORBI
         EIG = REIG(I)
         IF (EIG .LT. D0) THEN
            IF (PM%IPRLOC.GE.3 .AND. ABS(EIG).GT.VAR%THRZER) WRITE(LUPRI,9010) EIG
            EIG = D0
            REIG(I)   = D0
         END IF
         IF (EIG .LT. DTHRS) THEN

!           Use Taylor expansion to avoid division by zero and/or
!           to speed up evaluation (avoid call to sin/cos/sqrt).
!           For EIG < 1.0D-8 the 3. order Taylor expansion is correct
!             to machine precision.

            RCOS(1+IADDR) = D1-DP5*EIG
            RSIN(1+IADDR) = D1-DSIXTH*EIG
         ELSE
            SQEIG = SQRT(EIG)
            RCOS(1+IADDR) = COS(SQEIG)
            RSIN(1+IADDR) = SIN(SQEIG)/SQEIG
         END IF
      END DO
      IF (PM%IPRLOC.GE.10)THEN
         CALL HEADER('ROTCMO_MOD: sin(d)/d',-1)
         CALL PRQMAT(RSIN,NORBI,NORBI,NORBI,NORBI,NZ,IP,LUPRI)
         CALL HEADER('ROTCMO_MOD: cos(d)',-1)
         CALL PRQMAT(RCOS,NORBI,NORBI,NORBI,NORBI,NZ,IP,LUPRI)
      END IF
      CALL GETTIM(CPU3,WALL3)
 9010 FORMAT('WARNING from ROTCMO_MOD: Negative eigenvalue ',1P,D15.6)

!     Make V cos(d) V^{\dagger}
!     (use KTMP as temporatory storage)
!     ------------------------------------

      CALL QTRANS90_MOD('MOAO','S',D0,NORBI,NORBI,NORBI,NORBI, &
                        RTMP,NORBI,NORBI,NZ,IP,                &
                        RCOS,NORBI,NORBI,NZ,IP,                &
                        RV,NORBI,NORBI,NZ,IP,                  &
                        RV,NORBI,NORBI,NZ,IP,PM%IPRLOC)
      IF (PM%IPRLOC.GE.10)THEN
         CALL HEADER('ROTCMO_MOD: V cos(d) V^{dagger}',-1)
         CALL PRQMAT(RTMP,NORBI,NORBI,NORBI,NORBI,NZ,IP,LUPRI)
      END IF

!     Make V sin(d)/d V^{\dagger}
!     (use RCOS as temporatory storage)
!     ------------------------------------

      CALL QTRANS90_MOD('MOAO','S',D0,NORBI,NORBI,NORBI,NORBI, &
                        RCOS,NORBI,NORBI,NZ,IP,                &
                        RSIN,NORBI,NORBI,NZ,IP,                &
                        RV,NORBI,NORBI,NZ,IP,                  &
                        RV,NORBI,NORBI,NZ,IP,PM%IPRLOC)
      IF (PM%IPRLOC.GE.10)THEN
         CALL HEADER('ROTCMO_MOD: V sin(d)/d V^{dagger}',-1)
         CALL PRQMAT(RCOS,NORBI,NORBI,NORBI,NORBI,NZ,IP,LUPRI)
      END IF

!     Make V cos(d) V^{\dagger} - V sin(d)/d V^{\dagger} \kappa
!     ---------------------------------------------------------

      CALL QGEMM(NORBI,NORBI,NORBI,DM1,            &
                 'N','N',IP,RCOS,NORBI,NORBI,NZ,   &
                 'N','N',IP,XKAPPA,NORBT,NORBT,NZ, &
                 D1,IP,RTMP,NORBI,NORBI,NZ)
      IF (PM%IPRLOC.GE.10)THEN
         CALL HEADER('ROTCMO_MOD: exp(-kappa)',-1)
         CALL PRQMAT(RTMP,NORBI,NORBI,NORBI,NORBI,NZ,IP,LUPRI)
      END IF

!     Check if unitary

      CALL QGEMM(NORBI,NORBI,NORBI,D1,           &
                 'N','N',IP,RTMP,NORBI,NORBI,NZ, &
                 'H','N',IP,RTMP,NORBI,NORBI,NZ, &
                 D0,IP,RCOS,NORBI,NORBI,NZ)
      IOFF = 0
      I_MAX = 0
      J_MAX = 0
      D_MAX = d0
      DO I = 1,NORBI
      DO J = 1,NORBI
         IF (I.EQ.J) THEN
            IF (ABS(RCOS(1+IOFF)-D1).GT.D_MAX) THEN
               I_MAX = I
               J_MAX = J
               D_MAX = ABS(RCOS(1+IOFF)-D1)
            END IF
         ELSE
            IF (RCOS(1+IOFF).GT.D_MAX) THEN
               I_MAX = I
               J_MAX = J
               D_MAX = ABS(RCOS(1+IOFF))
            END IF
         END IF
         IOFF = IOFF + 1
      END DO
      END DO
      IF(PM%IPRLOC.GE.1)THEN
        IF(D_MAX.GT.VAR%THRZER)THEN
           WRITE(LUPRI,'(/A)')' WARNING: ROTCMO_MOD: maximum deviation from unitary matrix'
           WRITE(LUPRI,*)'WARNING: Element ',J_MAX,I_MAX,' value :',D_MAX
        ENDIF
      ENDIF

!     Calculate new orbitals: SCOEF_{old} exp(-kappa)
!     ---------------------------------------------

      IF (PM%IPRLOC.GE.10)THEN
         CALL HEADER('ROTCMO_MOD: old orbitals',-1)
         CALL PRQMAT(SCOEF,NBASI,NORBI,NBASI,NORBI,NZ,IP,LUPRI)
      END IF
!SK   CALL DZERO(RCOS,NORBXQ)
      CALL QGEMM(NBASI,NORBI,NORBI,D1,             &
                 'N','N',IP,SCOEF,NBASI,NORBI,NZ,  &
                 'N','N',IP,RTMP,NORBI,NORBI,NZ,   &
                 D0,IP,RCOS,NBASI,NORBI,NZ)
      CALL DCOPY(NBASI*NORBI*NZ,RCOS,1,SCOEF,1)
      IF (PM%IPRLOC.GE.5)THEN
         CALL HEADER('ROTCMO_MOD: new orbitals',-1)
         CALL PRQMAT(SCOEF,NBASI,NORBI,NBASI,NORBI,NZ,IP,LUPRI)
      END IF
      CALL GETTIM(CPU4,WALL4)
      IF (PM%IPRLOC.GE.2)THEN
         CPUTID = SECTID(CPU4-CPU1)
         WALLTID = SECTID(WALL4-WALL1)
         WRITE(LUPRI,9020) CPUTID,WALLTID
         CPUTID = SECTID(CPU3-CPU2)
         WALLTID = SECTID(WALL3-WALL2)
         WRITE(LUPRI,9021) CPUTID,WALLTID
      END IF
 9020 FORMAT(/'ROTCMO_MOD: total CPU (WALL) time       : ',A12,' (',A12,')')
 9021 FORMAT( 'ROTCMO_MOD: sin/sqrt/cos CPU (WALL) time: ',A12,' (',A12,')')

      DEALLOCATE( REIG )
      DEALLOCATE( RSIN )
      DEALLOCATE( RCOS )
      DEALLOCATE( RV   )
      DEALLOCATE( RTMP )

  END SUBROUTINE ROTCMO_MOD

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

      SUBROUTINE NUMLST_MOD(A,ILIST,NLIST,IRNGE1,IRNGE2,IFRP,IP)
!*****************************************************************************
!
!     Original purpose was to read the numbers in the string A which has
!     the format : '3..6, 8, -27..2, 123, 14..15'.
!     
!     !! Number 0 (Zero) NOT counted or stored !!
!
!     Input: 
!           IF(IP.NE.0) : Give list of numbers in array
!                         ILIST(NLIST).
!           ELSE : Only IP as OUTPUT
!
!     Output:
!           IP : number of elements in the list A.
!                Duplicates only counted once.
!
!     IRNGE1..IRNGE2 is range of allowed values in A.
!     ITMPAR   : Temporary storage
!     If the string starts with 'all', all elements in the range enters.
!
!                           Jon K. Laerdahl 30.12.96
!     Minor polish T.Saue Jan 7 1997
!
!     Additional feature (LV, 15-5-2001) : Determine list based on energy or
!     other criteria. This is more user-friendly for standard cases.
!
!     This option is toggled by starting the string by energy, followed by the
!     lower treshold, upper treshold and minimal gap (to prevent cutting through
!     quasi-degenerate orbitals).
!
!*****************************************************************************
      IMPLICIT NONE

      ! External variables
      INTEGER   :: IFRP, IP, NLIST, IRNGE1, IRNGE2
      CHARACTER :: A*(*)
      INTEGER   :: ILIST(NLIST)
  
      ! Internal variables
      LOGICAL :: TOBE
      INTEGER :: IER, NA, NORBT, LUCOEF, IOFF
      REAL*8  :: DUM

      REAL(8), ALLOCATABLE :: EIG(:)
      INTEGER, ALLOCATABLE :: IBG(:)
      INTEGER, ALLOCATABLE :: ITM(:)

  !------------------------------------------------------------------------------!
      NORBT = VAR%NORBT
      LUCOEF= VAR%LUCOEF

  !   Distinguish between the two branches

      NA = LEN(A)
      IF (NA.GT.80) NA = 80
      IF (INDEX(A(1:NA),'energy') .EQ. 0) THEN
         ALLOCATE( ITM(IRNGE2-IRNGE1+1), STAT=IER ); CALL ERRALL(IER,'NUMLST_MOD','ITM')
         CALL NUMLS1(A,ILIST,NLIST,IRNGE1,IRNGE2,IP,ITM)
         DEALLOCATE( ITM )
      ELSE
         ALLOCATE( EIG(NORBT), STAT=IER ); CALL ERRALL(IER,'NUMLST_MOD','EIG')
         ALLOCATE( IBG(NORBT), STAT=IER ); CALL ERRALL(IER,'NUMLST_MOD','IBG')
!
!        Get eigenvalues from DFCOEF
!
         INQUIRE(FILE='DFCOEF',EXIST=TOBE)
         IF(.NOT.TOBE.OR.VAR%I_DCBORB_SET.NE.1) THEN

!           DFCOEF does not (yet) exist, put IP = - 666, which makes 
!           NUMLS2 return with an empty list
!           Another showstopper occurs if we call this before the 
!           orbital common block is set (this info is currently needed
!           to correctly process the orbital energy list), this is
!           why we check for I_DCBORB_SET as well.


            IP = - 666
            CALL NUMLS2(A,ILIST,NLIST,IRNGE1,IRNGE2,IP,EIG)
         ELSE
            CALL OPNFIL(LUCOEF,'DFCOEF','OLD','PCMOUT')
            CALL REACMO(LUCOEF,'DFCOEF',DUM,EIG,IBG,DUM,4)
            CLOSE(LUCOEF,STATUS='KEEP')
            IOFF = VAR%IORB(IFRP)+VAR%NPSH(IFRP)
            CALL NUMLS2(A,ILIST,NLIST,IRNGE1,IRNGE2,IP,EIG(1+IOFF))
         ENDIF
         DEALLOCATE( EIG )
         DEALLOCATE( IBG )
      ENDIF

      END SUBROUTINE NUMLST_MOD

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine errall(ier,chsub,charr)
  implicit none

  integer      :: ier
  character(*) :: chsub, charr

  if(ier /= 0)then
    write(6,*)
    write(6,'(a,a,a,a)')' Error in allocation of ',charr,' in subroutine ',chsub
    write(6,'(a,i4)'   )' iostat =',ier
    stop 'ERROR in allocation (see output form more information)'
  endif

  end subroutine errall

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine get_overlap(smat)
  implicit none

  integer :: ier
  real*8  :: smat(:)

  ! open file
  open(60,file='OVERLAP.INP',status='old',form='formatted',access='sequential',iostat=ier)
  if(ier /= 0)then
    write(6,*) 'get_overlap: ERROR while opening OVERLAP.INP file'
    stop
  endif

  rewind 60

  read(60,*) smat

  ! close file
  close(60,status='keep')

  end subroutine get_overlap

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine selcfs_mod(cmo,ifrp,cbuf,lcb,jvec,npvec,nevec,knfbas,knorb)
  !==============================================================================!
  !          Pick out a set of vectors from CMO according to array JVEC          !
  !                     (based on Dirac routine SELCFS)                          !
  !==============================================================================!
  implicit none

  ! External variables
  integer, intent( in) :: ifrp, lcb, npvec, nevec, knfbas, knorb
  integer, intent( in) :: jvec(*)
  real*8,  intent( in) :: cmo(knfbas,knorb,var%nz)
  real*8,  intent(out) :: cbuf(knfbas,lcb,var%nz)

  ! Internal variables
  integer :: iz, nvec, i, ii

  !------------------------------------------------------------------------------!
  nvec = npvec + nevec
  do iz = 1,var%nz

    ! Positronic vectors
    do i = 1,npvec
      ii = var%npsh(ifrp)+1+jvec(i)
      call dcopy(var%nfbas(ifrp,0),cmo(1,ii,iz),1,cbuf(1,i,iz),1)
    enddo

    ! Electronic vectors
    do i = npvec+1,nvec
      ii = var%npsh(ifrp)+jvec(i)
      call dcopy(var%nfbas(ifrp,0),cmo(1,ii,iz),1,cbuf(1,i,iz),1)
    enddo
  enddo

  end subroutine selcfs_mod

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine qtrans90_mod(typ,trev,fadd,nrao,ncao,nrmo,ncmo, &
                          fao,lrao,lcao,nzao,iqao,           &
                          fmo,lrmo,lcmo,nzmo,iqmo,           &
                          tm1,lr1,lc1,nztm1,iqtm1,           &
                          tm2,lr2,lc2,nztm2,iqtm2,iprint)
  !==============================================================================!
  !       Like QTRANS, but uses fortran 90 memory allocation internally          !
  !      (based on Dirac QTRANS90 routine; uses standard f90 allocation)         !
  !==============================================================================!
  implicit none

  ! External variables
  character, intent(in) :: trev*1, typ*4
  integer,   intent(in) :: nrao, ncao, nrmo, ncmo, lrao, lcao, nzao, lrmo, lcmo, nzmo
  integer,   intent(in) :: lr1, lc1, nztm1, lr2, lc2, nztm2, iprint
  integer,   intent(in) :: iqao(*), iqmo(*), iqtm1(*), iqtm2(*)
  real*8,    intent(in) :: fadd, fao(*), fmo(*), tm1(*), tm2(*)

  ! Internal variables
  integer :: ier, nbuf, nzw, iqw(4)

  real*8, allocatable :: work(:)

  !------------------------------------------------------------------------------!
  !   find out how much buffer is needed in qtrans
  if(typ == 'AOMO')then
    if(ncao.ge.nrao)then
      call iqpack(iqao,nzao,iqtm2,nztm2,iqw,nzw)
      nbuf = nrao*ncmo*nzw
    else
      call iqpack(iqao,nzao,iqtm1,nztm1,iqw,nzw)
      nbuf = nrmo*ncao*nzw
    endif
  else
    if(ncmo >= nrmo)then
      call iqpack(iqmo,nzmo,iqtm2,nztm2,iqw,nzw)
      nbuf = nrmo*ncao*nzw
    else
      call iqpack(iqmo,nzmo,iqtm1,nztm1,iqw,nzw)
      nbuf = nrao*ncmo*nzw
    endif
  endif

  if(nbuf <= 0)then
   call quit(' *** error in QTRANS90: buffer space allocation not possible for length <= 0.***')
  end if

  allocate(work(nbuf), stat=ier); if(ier /= 0) call quit('qtrans90_mod: error in allocation !')
  work(:) = 0.0d0

  CALL qtrans(typ,trev,fadd,nrao,ncao,nrmo,ncmo, &
              fao,lrao,lcao,nzao,iqao,           &
              fmo,lrmo,lcmo,nzmo,iqmo,           &
              tm1,lr1,lc1,nztm1,iqtm1,           &
              tm2,lr2,lc2,nztm2,iqtm2,           &
              work,nbuf,iprint)

  deallocate(work, stat=ier); if (ier /= 0) call quit('qtrans90_mod: error in deallocation !')
      
  end subroutine qtrans90_mod

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine qdiag90_mod(nz,n,a,lra,lca,eig,matz,vec,lrv,lcv,ierr)
  !==============================================================================!
  !           Like QDIAG but uses f90 memory allocation internally               !
  !       (based on Dirac QDIAG90 routine; uses standard f90 allocation)         !
  !==============================================================================!
  implicit none

  ! External variables
  integer, intent(in) :: nz, n, lra, lca, matz, lrv, lcv, ierr
  real*8,  intent(in) :: a(*), eig(*), vec(*)

  ! Internal variables
  integer :: idum, iworkdum, info, llwork, lilwork
  real*8  :: ddum, workdum

  integer :: lenid, lwork, ier, m
  real*8, allocatable :: work(:)

  !------------------------------------------------------------------------------!

  lenid = 3 ! lenid used in MEMGET in QDIAG
  lwork = 4*lenid + 2*N
  if(nz > 1) lwork = lwork + N*NZ

#ifdef LAPACK_QDIAG
  if(matz == 1)then
    ddum = 0.0d0
    idum = 0
    workdum  = 0.0d0
    iworkdum = 0
    ! first let lapack estimate memory use
    info = 0
    call dsyevr('V','A','U',n,a,lra,ddum,ddum,idum,idum,0.0D0, &
                m,eig,vec,lrv,idum,workdum,-1,iworkdum,-1,info)
    llwork  = nint(workdum)
    lilwork = iworkdum
    lwork = llwork + lilwork + 2*n + 6*lenid
  endif
#endif

  if(lwork <= 0)then
   call quit(' *** error in qdiag90_mod: buffer space allocation not possible for length <= 0.***')
  endif

  allocate(work(lwork), stat=ier); if(ier /= 0) call quit('qdiag90_mod: error in allocation !')
  work(:) = 0.0d0

  call qdiag(nz,n,a,lra,lca,eig,matz,vec,lrv,lcv,work,lwork,ierr)

  deallocate(work, stat=ier); if (ier /= 0) call quit('qdiag90_mod: error in deallocation !')

  end subroutine qdiag90_mod

end module func_PipekMezey
