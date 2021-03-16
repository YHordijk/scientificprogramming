!      Copyright (c) 2018 by the authors of DIRAC.
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

  program vibcal
!*******************************************************************************!
!                                                                               !
!     For a given normal mode this program will calculate the                   !
!     expectation value for a property in the five lowest                       !
!     vibrational states                                                        !
!                                                                               !
!     Written by T. Saue Feb 13 2007                                            !
!     Further code by S. Komorovsky                                             !
!*******************************************************************************!
  use codata
  use num_der
  implicit none

#include "priunit.h"

  character :: potfil*100, prpfil*100, reply*1, strkey*30

  integer :: i, j, n, iprp, ipot, ilogg, lnblnk, ier, ierrpot, ierrprp, iang
  integer :: npoints, ipoint, itmp, nprint, dermeth, norder, ndim, iskip, keymeth

  integer, parameter :: nder=4

  real*8            :: u, v2, v3, v4, p0, p1, p2, um, omega, step
  real*8            :: fac, a0, a1, a2, a3, v, w, har, anh, tot, dif, dif2
  real*8            :: err1, err2, err3, edifb, edifh, edifb3
  real*8            :: chisq, yp, dev, xmin
  real*8            :: Rinf, Rinf3, Rinf4, dRinf, dRinf3, dRinf4, xmax
  real*8            :: total_00, total_01, total_02, total_03 

  real*8, parameter :: dp5=0.5d0, d1=1.0d0, dp16=0.0625d0, xtcm=2.19474625d+05
  real*8, parameter :: engerr=1.0d-10

  real*8, dimension(   0:4  ) :: b, h, b3, eb3, eb, eh
  real*8, dimension(nder  ,2) :: errder
  real*8, dimension(nder  ,3) :: derres
  real*8, dimension(2*nder,2) :: errwng

  real*8, allocatable, dimension(:) :: xcoor, ecoor, pcoor
  real*8, allocatable :: amat(:,:),bvec(:),cvec(:),dvec(:)

  real*8, external :: polval, polnder

  !-----------------------------------------------------------------------------!
  !                      Initialization of the code                             !
  !-----------------------------------------------------------------------------!
  errwng = 0.0d0  ! errwng(2*nder,2)
  ierrpot = 0; ierrprp = 0

  call titler('VIBCAL for properties','*',110)

  write(lupri,'(a)') 'Set print level [default:0]:'
  write(lupri,'(a)') '      0. Basic print level '
  write(lupri,'(a)') '      1. Print additional information'
  read(LUSTDIN,'(a)',iostat=ier) strkey; read(strkey,*,iostat=ier) nprint
  if(nprint /= 1) nprint = 0

  ! read reduced mass
  write(lupri,'(a)') 'Reduced mass (in Daltons)'
  read(LUSTDIN,*) u
  ! read the name of file with potential/property curve
  write(lupri,'(a)') 'Name of the input file (A30) with potential/property curve'
  read(LUSTDIN,'(a)') potfil
  write(lupri,'(a)')  'Select one of the following:'
  write(lupri,'(a)')  '  1. Bond lengths in Angstroms.'
  write(lupri,'(a)')  '  2. Bond lengths in atomic units.'
  write(lupri,'(a)')  ' (Note that all other quantities only in atomic units !)'
  read(LUSTDIN,*) iang
  write(lupri,'(a)')  'Select one of the following:'         
  write(lupri,'(a)')  '  1. Numerical derivatives (assumes fixed step length).'
  write(lupri,'(a)')  '  2. Polynomial fits             .'
  read(LUSTDIN,*) keymeth

  ! open input file
  ipot = 1
  open(ipot,file=potfil,status='old',form='formatted',access='sequential',iostat=ier)
  if(ier /= 0)then
    write(0,*) ' ERROR in opening input file POTFIL:',potfil
!gosia: this would write to unopened file, fixme
!    write(ilogg,*) ' ERROR in opening input file POTFIL:',potfil
    stop
  endif
  i=lnblnk(potfil)

  ! open output file
  ilogg = 2
  if (keymeth==1) then
    open(ilogg,file=potfil(1:i)//'.vibcal1', &
       status='unknown',form='formatted',access='sequential')
  elseif(keymeth==2)then  
    open(ilogg,file=potfil(1:i)//'.vibcal2', &
       status='unknown',form='formatted',access='sequential')
  else
    STOP 'Unknown option for derivative.'
  endif

  !-----------------------------------------------------------------------------!
  !                 Get force constants and property derivatives                !
  !-----------------------------------------------------------------------------!
  ! count number of points
  rewind ipot
  npoints = 0
  do
    read(ipot,*,iostat=ier)
    if(ier == -1)exit
    npoints = npoints + 1
  enddo
  if(nprint == 1)then
    write(ilogg,*)
    write(ilogg,'(a)') '  ******************************'
    write(ilogg,'(a)') '  ** DETAILED INFO ABOUT DER. **'
    write(ilogg,'(a)') '  ******************************'
    write(ilogg,'(a,i4,a)')'  We have ',npoints,' points'
  endif
  write(lupri,'(a,i4,a)')'We have ',npoints,' points'

  ! allocate arrays
  allocate ( xcoor(npoints),stat=ier ); if(ier /= 0 ) stop 'Problem with xcoor  all.'
  allocate ( ecoor(npoints),stat=ier ); if(ier /= 0 ) stop 'Problem with ecoor  all.'
  allocate ( pcoor(npoints),stat=ier ); if(ier /= 0 ) stop 'Problem with xcoor  all.'

  ! read points
  rewind ipot
  do i=1,npoints
    read(ipot,*,iostat=ier) xcoor(i), ecoor(i), pcoor(i)
    if(ier /= 0)then
      write(ilogg,*) ' ERROR during reading POTFIL:',potfil
      write(ilogg,*) ' STOP in vibcal'
      write(0,*) ' STOP in vibcal (see output)'
      stop
    endif
  enddo
  if(iang==1)then
    fac=d1/xtang
    call dscal(npoints,fac,xcoor,1)
  endif

  if (keymeth==1) then
    ! the derivation will be calculated in this point, with this step
    ipoint = (npoints / 2) + mod(npoints,2)
    step = xcoor(2) - xcoor(1)
    if(nprint==1)then
      write(ilogg,'(a,i4,a)')'  We use the mid-point (',ipoint,'th point) for calculation of derivatives'
      write(ilogg,'(a,es12.4)')'  step = ',step
      write(ilogg,*)
    endif

    ! calculate derivatives and estimate errors
    if(nprint == 1)write(ilogg,'(a)')' POTENTIAL DERIVATIONS'
    call nddrv(nder,'CENTRAL',ipoint,nprint,ilogg,npoints,xcoor,ecoor,engerr,derres,ierrpot)

    ! store results
    v2 = derres(2,1); v3 = derres(3,1); v4 = derres(4,1)
    xmax = xcoor(npoints)
    do i=1,nder
      errder(i,1) = max(derres(i,2),derres(i,3))
      errwng(i,1) = derres(i,2)
      errwng(i,2) = derres(i,3)
    enddo
    ! calculate derivatives and estimate errors
    if(nprint == 1)write(ilogg,'(a)')' PROPERTY DERIVATIONS'
    call nddrv(nder,'CENTRAL',ipoint,nprint,ilogg,npoints,xcoor,pcoor,0.0d0,derres,ierrprp)

    ! store results
    p0 = pcoor(ipoint); p1 = derres(1,1); p2 = derres(2,1)
    do i=1,nder
      errder(i,2) = max(derres(i,2),derres(i,3))
      errwng(i+nder,1) = derres(i,2)
      errwng(i+nder,2) = derres(i,3)
    enddo

  elseif(keymeth==2)then
    write(lupri,*) 'Give order of polynomial:'
    read(LUSTDIN,*) norder
    write(lupri,*) 'Select point for the calculation of derivatives'
    read(LUSTDIN,*)  xmin
    if(iang==1)then
      xmin=xmin/xtang
    endif
    do while (npoints.le.norder) 
      write(lupri,*) 'Too few points for a polynomial fit of order',norder
      write(lupri,*) '  Try again.'
      write(lupri,*) 'Give order of polynomial:'
      read(LUSTDIN,*) norder
    enddo
    ndim = norder + 1
  ! allocate arrays
    allocate ( amat(npoints,ndim),stat=ier ); if(ier /= 0 ) stop 'Problem with amat  all.'
    allocate ( bvec(npoints),stat=ier ); if(ier /= 0 ) stop 'Problem with bvec  all.'
    allocate ( cvec(ndim),stat=ier ); if(ier /= 0 ) stop 'Problem with cvec  all.'
    allocate ( dvec(ndim),stat=ier ); if(ier /= 0 ) stop 'Problem with dvec  all.'
  ! polynomial fit of potential energy curve
    call polsvd(ndim,npoints,amat,bvec,xcoor,ecoor,cvec,dvec,chisq,iskip)
    if(nprint == 1)then
      write(ilogg,'(72A1)') ('=',i=1,72)
      write(ilogg,'(3X,A,I3)') '* POTENTIAL CURVE: Polynomial fit of order:',norder
      write(ilogg,'(3X,A)') '* Coefficients:'
      write(ilogg,'(3X,A,I3,A,1P,E14.6)') ('c(',(I-1),'):  ',cvec(i),i=1,ndim)  
      write(ilogg,'(A)') '        X          Predicted E         Actual E            Relative error'
      do i = 1,npoints
        yp = polval(norder,cvec,xcoor(I))
        dev = (yp-ecoor(i))/ecoor(I)
        write(ilogg,'(3x,f6.3,4x,2e20.8,e16.4)') xcoor(I),yp,ecoor(i),dev
      enddo
      write(ilogg,'(3x,a,e9.4)') '* chi square :  ',chisq
      write(ilogg,'(3x,a,e9.4)') '* chi square per point:  ',chisq/npoints
      write(ilogg,'(3x,a,i3/)') '* Number of singularities(SVD): ',iskip
    endif
  ! calculate energy derivatives
    v2 = polnder(norder,cvec,xmin,2)
    v3 = polnder(norder,cvec,xmin,3)
    v4 = polnder(norder,cvec,xmin,4)
  ! polynomial fit of property curve
    call polsvd(ndim,npoints,amat,bvec,xcoor,pcoor,cvec,dvec,chisq,iskip)
    if(nprint == 1)then
      write(ilogg,'(72A1)') ('=',i=1,72)
      write(ilogg,'(3X,A,I3)') '* PROPERTY CURVE: Polynomial fit of order:',norder
      write(ilogg,'(3X,A)') '* Coefficients:'
      write(ilogg,'(3X,A,I3,A,1P,E14.6)') ('c(',(I-1),'):  ',cvec(i),i=1,ndim)  
      write(ilogg,'(A)') '        X          Predicted P         Actual P            Relative error'
      do i = 1,npoints
        yp = polval(norder,cvec,xcoor(I))
        dev = (yp-pcoor(i))/ecoor(I)
        write(ilogg,'(3x,f6.3,4x,2e20.8,e16.4)') xcoor(I),yp,pcoor(i),dev
      enddo
      write(ilogg,'(3x,a,e9.4)') '* chi square :  ',chisq
      write(ilogg,'(3x,a,e9.4)') '* chi square per point:  ',chisq/npoints
      write(ilogg,'(3x,a,i3/)') '* Number of singularities(SVD): ',iskip
    endif
  ! calculate property derivatives
    p0 = polval(norder,cvec,xmin)
    p1 = polnder(norder,cvec,xmin,1)
    p2 = polnder(norder,cvec,xmin,2)
  ! deallocate
      deallocate ( amat )
      deallocate ( bvec )
      deallocate ( cvec )
      deallocate ( dvec )
  else
    STOP 'Unknown option for derivative.'
  endif

  write(ilogg,*)
  write(ilogg,'(a)') '  ******************************'
  write(ilogg,'(a)') '  ****  SOME INFORMATION    ****'
  write(ilogg,'(a)') '  ******************************'

  write(ilogg,*) 'Reduced mass in Daltons: ',u
  write(ilogg,*) '2nd potential derivative: ', v2
  write(ilogg,*) '3rd potential derivative: ', v3
  write(ilogg,*) '4th potential derivative: ', v4

  write(ilogg,*) '0th property derivative: ', p0
  write(ilogg,*) '1st property derivative: ', p1
  write(ilogg,*) '2nd property derivative: ', p2

  !-----------------------------------------------------------------------------!
  !                          Calculate PV shift                                 !
  !-----------------------------------------------------------------------------!
  !--------------------------------------------!
  ! calculate harmonic and anharmonic PV shift !
  !--------------------------------------------!
  um = u*xfamu
  omega = sqrt(v2/um)
  fac=d1/um/omega
  a0 = p0
  a1 = dp5*p2*fac
  a2 = -dp5*p1*v3*fac*fac/omega
  a3 = -dp16*p2*v4*fac*fac*fac/omega

  if(nprint == 1)then
    write(ilogg,*)
    write(ilogg,'(a)') '  ******************************'
    write(ilogg,'(a)') '  ***  CONTRIBUTIONS TO PV   ***'
    write(ilogg,'(a)') '  ******************************'
    write(ilogg,'(a,a)') '    n  Term 1 (P0)     Term 2 (P2)     Term 3 (P1 V3) ', &
                         ' Term 4 (P2 V4)  Harmonic cont.  Anharm. cont.   Total value'
  endif
  do n = 0,4
    v = n + dp5
    w = n*n+n+dp5
    har = a0 + a1*v
    anh = a2*v + a3*w
    tot = har + anh 
    b(n) = tot         
    b3(n) = har + a2*v
    h(n) = har
    if(nprint == 1) write(ilogg,'(i5,7e16.8)') n,a0,a1*v,a2*v,a3*w,har,anh,tot
  enddo

  !------------------!
  ! error estimation !
  !------------------!
  if(keymeth == 1)then
    err1 = dp5*fac*dabs(errder(2,2))
    err2 = dp5*fac*fac/omega * (dabs(p1*errder(3,1)) + dabs(v3*errder(1,2)))
    err3 = dp16*fac*fac*fac/omega * (dabs(p2*errder(4,1)) + dabs(v4*errder(2,2)))
  endif

  if(keymeth == 1)then
    if(nprint == 1)then
      write(ilogg,*)
      write(ilogg,'(a)') '  ******************************'
      write(ilogg,'(a)') '  *  CONTRIBUTIONS TO ERROR    *'
      write(ilogg,'(a)') '  ******************************'
      write(ilogg,'(a,a)') '    n  Term 1 (P0)     Term 2 (P2)     Term 3 (P1 V3) ', &
                           ' Term 4 (P2 V4)  Total error'
    endif
    do n = 0,4
      v = n + dp5
      w = n*n+n+dp5
      eb(n) = v*err1 + v*err2 + w*err3
      eb3(n) = v*err1 + v*err2 
      eh(n) = v*err1
      if(nprint == 1) write(ilogg,'(i5,5e16.8)') n,0.0d0,v*err1,v*err2,w*err3,eb(n)
    enddo
    edifb  = xthz * 1000.0d0 * (eb(0) + eb(1))
    edifb3 = xthz * 1000.0d0 * (eb3(0) + eb3(1))
    edifh  = xthz * 1000.0d0 * (eh(0) + eh(1))
  endif

  !---------------!
  ! Print results !
  !---------------!
  write(ilogg,*)
  write(ilogg,'(a)') '  ******************************'
  write(ilogg,'(a)') '  ********   RESULTS    ********'
  write(ilogg,'(a)') '  ******************************'

  dif = b(1)-b(0)
  dif2 = dif*xthz*1000.0d0
  if(keymeth == 1)then
    write(ilogg,'(a,e40.8,a,f40.1,a,f40.4)')'  Difference :',dif, ' au = ',dif2, ' mHz  +-',edifb
    write(ilogg,'(a,f40.3,a,f40.4)')'  PV shift:',2.0d0*dif2,' mHz  +-',2.0d0*edifb
  else
    write(ilogg,'(a,e40.8,a,f40.1,a)')'  Difference :',dif, ' au = ',dif2, ' mHz'
    write(ilogg,'(a,f40.3,a)')'  PV shift:',2.0d0*dif2,' mHz'
  endif

  if(keymeth == 1)then
    dif = b3(1)-b3(0)
    dif2 = dif*xthz*1000.0d0
    write(ilogg,'(a)')'  Cubic model:'
    write(ilogg,'(a,e40.8,a,f40.1,a,f40.4)')'  Difference :',dif, ' au = ',dif2, ' mHz  +-',edifb3
    write(ilogg,'(a,f40.3,a,f40.4)')'  PV shift:',2.0d0*dif2,' mHz  +-',2.0d0*edifb3
  endif

  write(ilogg,'(a)')'  Harmonic model:'
  dif = h(1)-h(0)
  dif2 = dif*xthz*1000.0d0
  if(keymeth == 1)then
    write(ilogg,'(a,e40.8,a,f40.1,a,f40.4)')'  Difference :',DIF, ' au = ',dif2, ' mHz  +-',edifh
    write(ilogg,'(a,f40.3,a,f40.4)')'  PV shift:',2.0d0*dif2,' mHz  +-',2.0d0*edifh
  else
    write(ilogg,'(a,e40.8,a,f40.1,a)')'  Difference :',DIF, ' au = ',dif2, ' mHz'
    write(ilogg,'(a,f40.3,a)')'  PV shift:',2.0d0*dif2,' mHz'
  endif

  !-----------------------------------------------------------------------------!
  !              Hint if the choice of points were correct                      !
  !        We assume that the potential curve is Morse potential                !
  !             Rinf is then inflex point of Morse potential                    !
  !-----------------------------------------------------------------------------!
  if(keymeth == 1)then
    Rinf3  = dlog(2.0d0) * 3.0d0 * dabs(v2/v3)
    dRinf3 = dlog(2.0d0) * 3.0d0 * (dabs(v3*errder(2,1)) + dabs(v2*errder(3,1))) / (v3 * v3)

    Rinf4  = dlog(2.0d0) * dsqrt(7.0d0) * dsqrt(v2/v4)
    dRinf4 = (dabs(v4*errder(2,1)) + dabs(v2*errder(4,1))) / (v4 * v4)
    dRinf4 = dlog(2.0d0) * dsqrt(7.0d0) * dRinf4 / (2.0d0 * dsqrt(v2/v4))

    Rinf  = Rinf3;  if(dRinf4 < dRinf3) Rinf  = Rinf4
    dRinf = dRinf3; if(dRinf4 < dRinf3) dRinf = dRinf4

    !-----------------!
    ! WARNING section !
    !-----------------!
    write(ilogg,*)
    write(ilogg,'(a)') '  ******************************'
    write(ilogg,'(a)') '  ********   WARNINGS   ********'
    write(ilogg,'(a)') '  ******************************'

    ! compare different errors
    write(ilogg,'(a,a)') '# Compare errors comming from: error in f(x) (err_func)', &
                         ' and truncating of expansion (err_trun)'
    itmp = 0
    do i=1, 2*nder
      if(errwng(i,2) > errwng(i,1))then
        if(itmp==0)write(ilogg,'(6x,a)') 'WARNING: The err_func is greater than err_trun '
        itmp = 1
        if(i<nder)then
          write(ilogg,'(15x,a,i2,a,d10.3,a,d10.3)') &
            'Potential',i,'. derivation:  err_func =',errwng(i,2),'   err_trun =',errwng(i,1)
        else
          write(ilogg,'(15x,a,i2,a,d10.3,a,d10.3)') &
            'Property',i,'. derivation:  err_func =',errwng(i,2),'   err_trun =',errwng(i,1)
        endif
      endif
    enddo
    if(itmp==1)write(ilogg,'(6x,a)')'CHECK convergence of numerical derivations!!!'
    if(itmp==0)write(ilogg,'(a)')'  ************  OK  ************'; write(ilogg,*)

    ! check inflex point
    write(ilogg,'(a,a)') '# Check the estimated inflex point Rinf (based on Morse potential)'
    itmp = 0
    if(Rinf < xmax)then
      itmp = 1
      write(ilogg,'(6x,a)') 'WARNING: Your interval <-I,I> is TOO big: I > Rinf '
      write(ilogg,'(6x,a)') '         Your derivations have probably bad convergence'
    endif
    if((itmp==0).and.((Rinf/2.0d0) < xmax))then
      itmp = 1
      write(ilogg,'(6x,a)') 'WARNING: Your interval <-I,I> is probably too big: I > Rinf/2 '
    endif
    if(itmp==1)then
      write(ilogg,'(15x,a,d12.4)') 'I    =',xmax
      write(ilogg,'(15x,a,d12.4,a,d12.4)') 'Rinf =',Rinf,'  +-',dRinf
      do i=1,npoints
        if( ((Rinf-dRinf)/2.0d0) < xcoor(i) )then
          write(ilogg,'(6x,a,a,i3,a,d11.4)')'TRY to recalculate PV on interval',&
            ' <-x_n,...,x_0,...,x_n>:  n =',i-1-ipoint,';  x_n =',xcoor(i-1)
          exit
        endif
      enddo
    endif
    if(itmp==1)write(ilogg,'(6x,a)')'CHECK convergence of numerical derivations!!!'
    if(itmp==0)write(ilogg,'(a)')'  ************  OK  ************'; write(ilogg,*)

    ! check monotonousness
    write(ilogg,'(a,a)')'# Check if the convergence of derivations are monotonous'
    itmp = 0
    if(ierrpot /= 0)then
      itmp = 1
      write(ilogg,'(6x,a)')'WARNING: The convergence of potential der. is not monotonous '
      write(ilogg,'(15x,a,i5,a)')'Check these derivations:',ierrpot,' ("4321", 1-true, 0-false)'
    endif
    if(ierrpot /= 0)then
      itmp = 1
      write(ilogg,'(6x,a)')'WARNING: The convergence of property der. is not monotonous '
      write(ilogg,'(15x,a,i5,a)')'Check these derivations:',ierrpot,' ("4321", 1-true, 0-false)'
    endif
    if(itmp==1)write(ilogg,'(6x,a)')'CHECK convergence of numerical derivations!!!'
    if(itmp==0)write(ilogg,'(a)')'  ************  OK  ************'; write(ilogg,*)

  endif

  i=lnblnk(potfil)
  write(lupri,*) 'Results have been written to ',potfil(1:i)//'.vibcal'

  !-----------------------------------------------------------------------------!
  deallocate ( xcoor  )
  deallocate ( ecoor  )
  deallocate ( pcoor  )

  close(ilogg, status='keep')
  close(ipot, status='keep')
  !-----------------------------------------------------------------------------!
  end program
