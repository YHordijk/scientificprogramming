module num_der
implicit none

public  :: nddrv

private :: forwder, backder, centder, mul2exp, stir1, forwn, backn, centn, centnmu

integer, private :: luout

contains

!================================================================================!
!================================================================================!

  subroutine nddrv(nder,derstr,ipoint,nprint,lunit,npoints,xcoor,ycoor,yerr,derres,ierr)
  !------------------------------------------------------------------------------!
  !      Driver for numerical calculation of derivatives in some point           !
  !------------------------------------------------------------------------------!
  !  nder    -->  inp: maximum order of calculated derivation                    !
  !  derstr  -->  inp: what kind of method you want to use                       !
  !          -->  FORWARD   -->  eq. 5.3.7 in [1]                                !
  !          -->  BACKWARD  -->  eq. 5.3.9 in [1]                                !
  !          -->  CENTRAL   -->  eq. 5.3.10 - 5.3.16 in [1]                      !
  !  ipoint  -->  inp: pointer on point in which will be calc. the derivation    !
  !  nprint  -->  inp: print level                                               !
  !          -->  0  -->  do not print results                                   !
  !          -->  1  -->  print results                                          !
  !  xcoor(npoints) -->  inp: x coordinates, they must have uniform step length  !
  !  ycoor(npoints) -->  inp: y=f(x); functional values in points xcoor          !
  !  yerr           -->  inp: error in calculation of ycoor                      !
  !  derres(nder,1) -->  out: array containing all derivations                   !
  !  derres(nder,2) -->  out: error estimation from truncating of expansion      !
  !  derres(nder,3) -->  out: error estimation from error in ycoor               !
  !------------------------------------------------------------------------------!
  !     [1]:  F.B. Hildebrand, Introduction to Numerical Analysis, DPI, Inc.,    !
  !           New York, 1974.                                                    !
  !------------------------------------------------------------------------------!
! use unit_test_generator
  implicit none

  ! external variables
  character(*), intent(in) :: derstr
  integer,      intent(in) :: nder, nprint, lunit, npoints, ipoint
  integer,      intent(out):: ierr
  real(8),      intent(in) :: yerr

  real(8),      intent(in), dimension(npoints) :: xcoor, ycoor
  real(8),      intent(out),dimension(nder, 3) :: derres

  ! internal variables
  character :: strres*14

  integer :: i, k, ntmp, itmp, ier

  real(8), parameter :: tresh = 1.0d-10
  real(8) :: step, tmp

  integer, allocatable, dimension(:)   :: ndmax
  real(8), allocatable, dimension(:,:) :: dcon, err3c, dpar
  integer, save :: icount = 0

  !------------------------------------------------------------------------------!
  !           Write to the file input information for unit testing               !
  !------------------------------------------------------------------------------!
#ifdef CREATE_UNIT_TEST_nddrv
  icount = icount + 1
  call set_suffix(icount)
  call set_driver_routine('nddrv')
  call add_module('num_der')
 
  call add_inp_variable(nder   ,'nder   ')
  call add_inp_variable(derstr ,'derstr ')
  call add_inp_variable(ipoint ,'ipoint ')
  call add_inp_variable(nprint ,'nprint ')
  call add_inp_variable(lunit  ,'lunit  ')
  call add_inp_variable(npoints,'npoints')
  call add_inp_variable(xcoor  ,'xcoor  ')
  call add_inp_variable(ycoor  ,'ycoor  ')
  call add_inp_variable(yerr   ,'yerr   ')
  call add_inp_variable(ierr   ,'ierr   ')
#endif

  !------------------------------------------------------------------------------!
  luout = lunit
  step = dabs(xcoor(1)-xcoor(2))

  ! check the uniformity of step length
  do i=2, npoints-1
    tmp = dabs(xcoor(i)-xcoor(i+1))
    if(dabs(step-tmp) > tresh)then
      write(luout,*)" You don't have uniform step length"
      write(luout,'(a,d20.12)')'  x(1)-x(2) =',step
      write(luout,'(a,i3,a,d20.12)')'  i =',i,';  x(i)-x(i+1) =',tmp
      write(luout,*)' STOP in module num_der'
      write(0,*)' STOP in module num_der (check output)'
      stop
    endif
  enddo

  ! allocate arrays
  allocate ( ndmax(    nder    ),stat=ier ); if(ier /= 0 ) stop 'Problem with ndmax  all.'
  allocate ( dcon (nder,npoints),stat=ier ); if(ier /= 0 ) stop 'Problem with dcon   all.'
  allocate ( err3c(nder,npoints),stat=ier ); if(ier /= 0 ) stop 'Problem with err3c  all.'
  allocate ( dpar (nder,npoints),stat=ier ); if(ier /= 0 ) stop 'Problem with dpar   all.'

  ndmax  = 0     ! ndmax (nder)
  derres = 0.0d0 ! derres(nder,3)

  !------------------------------------------------------------------------------!
  !                           Calculate derivations                              !
  !------------------------------------------------------------------------------!
  ! calculate derivatives using forward differences
  if(derstr == 'FORWARD')then
    do i=1,nder
      call forwder(i,ipoint,step,ycoor,dcon(i,:),err3c(i,:),ndmax(i))
    enddo
  endif

  ! calculate derivatives using backward differences
  if(derstr == 'BACKWARD')then
    do i=1,nder
      call backder(i,ipoint,step,ycoor,dcon(i,:),err3c(i,:),ndmax(i))
    enddo
  endif

  ! calculate derivatives using central differences
  if(derstr == 'CENTRAL')then
    do i=1,nder
      call centder(i,ipoint,step,ycoor,dcon(i,:),err3c(i,:),ndmax(i))
    enddo
  endif

  ! calculate result from contributions
  do i=1, nder
    itmp=0
    if(ndmax(i) /= 0)then
      do k=1,ndmax(i)
        derres(i,1) = derres(i,1) + dcon(i,k)
        dpar(i,k)   = derres(i,1)
        ntmp        = int(dlog10(dabs(dcon(i,k))))  ! get mantissa of dcon(i,k)
        err3c(i,k)  = err3c(i,k) * yerr
        derres(i,3) = derres(i,3) + err3c(i,k)
        if(k/=1)then
          if(dabs(dcon(i,k)) > dabs(dcon(i,k-1)))itmp=1
        endif
      enddo
      if(itmp==1)ierr = ierr + 10**(i-1)
      derres(i,2) = dexp(dfloat(ntmp)*dlog(10.0d0))
    endif
  enddo

  ! Print results
  if(nprint == 1)then
    write(luout,'(25x,a,10x,a,5x,a)') 'Derivation','Partial sum','Err. from err. in f(x)'
    do i=1, nder
      if(derstr == 'CENTRAL') strres = ' SUM delta(mu)'
      if(derstr == 'FORWARD') strres = ' SUM forward  '
      if(derstr == 'BACKWARD')strres = ' SUM backward '
      if(ndmax(i) /= 0)then
        do k=1,ndmax(i)
          if(k==1)then
            write(luout,'(1x,i1,a,i1,a,es21.10,21x,es21.10)') &
             k,'. order d^',i,'  = ',dcon(i,k),err3c(i,k)
          else
            write(luout,'(1x,i1,a,i1,a,3es21.10)') &
              k,'. order d^',i,'  = ',dcon(i,k),dpar(i,k-1),err3c(i,k)
          endif
        enddo
        write(luout,'(a14,a,21x,2es21.10)')strres,' = ',derres(i,1),derres(i,3)
        write(luout,*)
      endif
    enddo
  endif

  ! deallocate arrays
  deallocate( ndmax )
  deallocate( dcon  )
  deallocate( err3c )
  deallocate( dpar  )

  !------------------------------------------------------------------------------!
  !           Write to the file ouput information for unit testing               !
  !------------------------------------------------------------------------------!
#ifdef CREATE_UNIT_TEST_nddrv
  call add_out_variable(derres,'derres')
  call generate_test()
#endif

  !------------------------------------------------------------------------------!
  end subroutine nddrv

!================================================================================!
!================================================================================!
  subroutine forwder(ider,ipoint,space,ycoor,dcon,derrc,ndmax)
  !------------------------------------------------------------------------------!
  !           Calculate numerical derivatives using forward differences          !
  !                       We use the formula 5.3.7 in [1]                        !
  !------------------------------------------------------------------------------!
  ! ider   --> inp: order of derivation                                          !
  ! ipoint --> inp: pointer on point in which will be calculated the derivation  !
  ! space  --> inp: step beteen points                                           !
  ! ycoor  --> inp: functional values of points                                  !
  ! dcon   --> out: individual contributions of the derivation                   !
  ! derrc  --> out: if the points (ycoor) are calculated with some error "err",  !
  !                 the estimated error of individual contributions is derrc*err !
  ! ndmax  --> out: number of individual contributions                           !
  !------------------------------------------------------------------------------!
  implicit none

  ! external variables
  integer, intent(in) :: ider, ipoint
  integer             :: ndmax
  real(8), intent(in) :: space

  real(8), dimension(:)  :: dcon, derrc, ycoor

  ! internal variables
  integer :: i, j, ier, npoints, nmax
  real(8) :: tmp1, tmp2, xstep

  integer, allocatable, dimension(:,:) :: istirc1
  real(8), allocatable, dimension(:)   :: tp
  !------------------------------------------------------------------------------!
  !                             Initialize code                                  !
  !------------------------------------------------------------------------------!
  !--------------------!
  !     Constants      !
  !--------------------!
  npoints = ubound(ycoor,1)  ! get the dimension of ycoor array

  nmax = npoints - ipoint    ! maximum number of differences
  if(nmax < ider)then
    write(luout,'(a,i3,a)') "You don't have enough points for calculating ",ider,"th derivation"
    write(luout,*) 'WARNING in forwder (module num_der)'
    ndmax = 0
    return
  endif

  ! step factor
  xstep = 1.0d0
  do i=1,ider
    xstep = xstep * space
  enddo

  !----------------------!
  !   allocate fields    !
  !----------------------!
  allocate ( istirc1(0:nmax,0:ider),stat=ier ); if(ier /= 0 ) stop 'Problem with istirc1 all.'
  allocate ( tp     (    nmax     ),stat=ier ); if(ier /= 0 ) stop 'Problem with tp      all.'

  ! set array tp
  tp = 1.0d0  ! tp(ntp)
  do i=1,nmax
    do j=1,i
      tp(i) = tp(i) * 2.0d0
    enddo
!   tp(i) = dsqrt(tp(i))
  enddo

  ! Stirling numbers of the first kind
  call stir1(nmax,ider,istirc1)

  !------------------------------------------------------------------------------!
  !                                 Main code                                    !
  !------------------------------------------------------------------------------!
  tmp1 = 1.0d0
  do i=1, nmax-ider+1
    tmp2 = dfloat(istirc1(ider+i-1,ider)) / tmp1
    if(mod(i-1,2) == 1)tmp2 = - tmp2
    dcon(i) = tmp2 * forwn(ider+i-1,ycoor,ipoint,npoints,ier) / xstep
    tmp1 = tmp1 * dfloat(ider + i)

    derrc(i) = tp(ider+i-1) * dabs(tmp2) / xstep
  enddo
  ndmax = nmax-ider+1

  deallocate( istirc1 )
  deallocate( tp      )

  !------------------------------------------------------------------------------!
  end subroutine forwder

  !==============================================================================!
  !==============================================================================!

  subroutine backder(ider,ipoint,space,ycoor,dcon,derrc,ndmax)
  !------------------------------------------------------------------------------!
  !          Calculate numerical derivatives using backward differences          !
  !                       We use the formula 5.3.9 in [1]                        !
  !------------------------------------------------------------------------------!
  !           For description of variables see subroutine forwder                !
  !------------------------------------------------------------------------------!
  implicit none

  ! external variables
  integer, intent(in) :: ider, ipoint
  integer             :: ndmax
  real(8), intent(in) :: space

  real(8), dimension(:)  :: dcon, derrc, ycoor

  ! internal variables
  integer :: i, j, ier, npoints, nmax
  real(8) :: tmp1, tmp2, xstep

  integer, allocatable, dimension(:,:) :: istirc1
  real(8), allocatable, dimension(:)   :: tp
  !------------------------------------------------------------------------------!
  !                             Initialize code                                  !
  !------------------------------------------------------------------------------!
  !--------------------!
  !     Constants      !
  !--------------------!
  npoints = ubound(ycoor,1)  ! get the dimension of ycoor array

  nmax = ipoint - 1          ! maximum number of differences
  if(nmax < ider)then
    write(luout,'(a,i3,a)') "You don't have enough points for calculating ",ider,"th derivation"
    write(luout,*) 'WARNING in backder (module num_der)'
    ndmax = 0
    return
  endif

  xstep = 1.0d0
  do i=1,ider
    xstep = xstep * space
  enddo

  !----------------------!
  !   allocate fields    !
  !----------------------!
  allocate ( istirc1(0:nmax,0:ider),stat=ier ); if(ier /= 0 ) stop 'Problem with istirc1 all.'
  allocate ( tp     (    nmax     ),stat=ier ); if(ier /= 0 ) stop 'Problem with tp      all.'

  ! set array tp
  tp = 1.0d0  ! tp(ntp)
  do i=1,nmax
    do j=1,i
      tp(i) = tp(i) * 2.0d0
    enddo
!   tp(i) = dsqrt(tp(i))
  enddo

  ! Stirling numbers of the first kind
  call stir1(nmax,ider,istirc1)

  !------------------------------------------------------------------------------!
  !                                 Main code                                    !
  !------------------------------------------------------------------------------!
  tmp1 = 1.0d0
  do i=1, nmax-ider+1
    tmp2 = dfloat(istirc1(ider+i-1,ider)) / tmp1
    dcon(i) = tmp2 * backn(ider+i-1,ycoor,ipoint,npoints,ier) / xstep
    tmp1 = tmp1 * dfloat(ider + i)

    derrc(i) = tp(ider+i-1) * dabs(tmp2) / xstep
  enddo
  ndmax = nmax-ider+1

  deallocate( istirc1 )
  deallocate( tp      )

  !------------------------------------------------------------------------------!
  end subroutine backder

  !==============================================================================!
  !==============================================================================!

  subroutine centder(ider,ipoint,space,ycoor,dcon,derrc,ndmax)
  !------------------------------------------------------------------------------!
  !       Calculate even and odd derivations using central diferences            !
  !                We use the formulas 5.3.10 - 5.3.16 in [1]                    !
  !------------------------------------------------------------------------------!
  !           For description of variables see subroutine forwder                !
  !------------------------------------------------------------------------------!
  implicit none

  ! external variables
  integer, intent(in) :: ider, ipoint
  integer             :: ndmax
  real(8), intent(in) :: space

  real(8), dimension(:) :: dcon, derrc, ycoor

  ! internal variables
  integer :: npoints, ntp, nmax, npoints2, ipoint2
  integer :: i, ii, j, ier
  real(8) :: tmp1, tmp2, xstep

  real(8), allocatable, dimension(:) :: expan, expan1, expans, ares, ycoor2, tp

  !------------------------------------------------------------------------------!
  !                             Initialize code                                  !
  !------------------------------------------------------------------------------!
  !--------------------!
  !     Constants      !
  !--------------------!
  npoints = ubound(ycoor,1)  ! get the dimension of ycoor array

  nmax = 2 * min(npoints-ipoint,ipoint-1)
  if(nmax < ider)then
    write(luout,'(a,i3,a)') "You don't have enough points for calculating ",ider,"th derivation"
    write(luout,*) 'WARNING in centder (module num_der)'
    ndmax = 0
    return
  endif

  ! get dimension for tp array
  if(mod(ider,2) == 0) ntp = nmax-ider+2
  if(mod(ider,2) == 1) ntp = nmax-ider+1

  ! step of central differences is half the distance between points
  ! therefore we reorganize n points (ycoor) to 2n-1 points (ycoor2)
  npoints2 = npoints * 2 - 1
  ipoint2 = ipoint * 2 - 1

  ! step factor
  xstep = 1.0d0
  do i=1,ider
    xstep = xstep * space
  enddo

  !----------------------!
  !   allocate fields    !
  !----------------------!
  allocate ( expan (0:nmax)  ,stat=ier ); if(ier /= 0 ) stop 'Problem with expan  all.'
  allocate ( expan1(0:nmax)  ,stat=ier ); if(ier /= 0 ) stop 'Problem with expan1 all.'
  allocate ( expans(0:nmax)  ,stat=ier ); if(ier /= 0 ) stop 'Problem with expans all.'
  allocate ( ares  (0:nmax)  ,stat=ier ); if(ier /= 0 ) stop 'Problem with ares   all.'
  allocate ( ycoor2(npoints2),stat=ier ); if(ier /= 0 ) stop 'Problem with ycoor2 all.'
  allocate ( tp    (ntp     ),stat=ier ); if(ier /= 0 ) stop 'Problem with tp     all.'

  !  insert zero between points
  ii = 0
  do i=1, npoints2
    if(mod(i,2)==0)then
      ycoor2(i) = 0.0d0
    else
      ii = ii + 1
      ycoor2(i) = ycoor(ii)
    endif
  enddo

  ! set array tp
  tp = 1.0d0  ! tp(ntp)
  do i=1,ntp
    do j=1,i
      tp(i) = tp(i) * 2.0d0
    enddo
!   tp(i) = dsqrt(tp(i))
  enddo

  !------------------------------------------------------------------------------!
  !                                 Main code                                    !
  !------------------------------------------------------------------------------!
  ! Expansion coefficients of Inverse Hyperbolic Sine
  ! Every term is aditionaly multiplied with 2^(n-1)
  ! f(x) = 2 * sinh^(-1)(x/2) 
  !      = sum_{n=0} { (-1)^n  [(2n-1)!!]^2 / [ (2n+1)! 2^{2n} ] * x^{2n+1} }
  tmp1 = 1.0d0
  tmp2 = 1.0d0
  expan(0) = 0.0d0; expans(0) = 0.0d0
  expan(1) = 1.0d0; expans(1) = 1.0d0
  do i=2, nmax
    tmp2 = tmp2 * dfloat(i) * 2.0d0
    if(mod(i,2) == 0)then
      expan(i) = 0.0d0
      expans(i) = 0.0d0
    else
      tmp1 = - tmp1 * dfloat((i-2)*(i-2))
      expan(i) = tmp1 / tmp2
      expans(i) = tmp1 / tmp2
    endif
  enddo

  !---------------------------------------------------!
  !  Get expansion coefficients for even derivatives  !
  !---------------------------------------------------!
  do i=1, ider/2
    call mul2exp(nmax,nmax,nmax,expan,expan,ares)
    do j=0, nmax
      expan(j) = ares(j)
    enddo
  enddo

  !--------------------------------------------------!
  !  Get expansion coefficients for odd derivatives  !
  !--------------------------------------------------!
  if(mod(ider,2) == 1)then
    ! Expansion coefficients of function:
    ! f(x) = 1 / sqrt(1+x*x/4)
    !      = sum_{n=0} { (-1)^n  [(2n-1)!!]^2 / [ (2n)! 2^{2n} ] * x^{2n} }
    tmp1 = 1.0d0
    tmp2 = 2.0d0
    expan1(0) = 1.0d0
    expan1(1) = 0.0d0
    do i=2, nmax
      tmp2 = tmp2 * dfloat(i) * 2.0d0
      if(mod(i,2) == 1)then
        expan1(i) = 0.0d0
      else
        tmp1 = - tmp1 * dfloat((i-1)*(i-1))
        expan1(i) = tmp1 / tmp2
      endif
    enddo

    ! Expansion coefficients of function:
    ! f(x) = ( 1 / sqrt(1+x*x/4) ) * 2 * sinh^(-1)(x/2)
    call mul2exp(nmax,nmax,nmax,expan1,expans,ares)

    ! Get final expansion coefficients
    if(ider /= 1)then
      do j=0, nmax
        expan1(j) = ares(j)
      enddo
      call mul2exp(nmax,nmax,nmax,expan1,expan,ares)
    endif

  endif

  !---------------------------------------------------------!
  !   Calculate individual terms and estimate their errors  !
  !---------------------------------------------------------!
  do i=1,(nmax-ider)/2 + 1
    if(mod(ider,2) == 0)then
      dcon(i)  = ares(2*i) * centn(2*i,ycoor2,ipoint2,npoints2,ier) / xstep
      derrc(i) = tp(2*i) * dabs(ares(2*i)) / xstep
    else
      dcon(i)  = ares(2*i-1) * centnmu(2*i-1,ycoor2,ipoint2,npoints2,ier) / xstep
      derrc(i) = tp(2*i-1) * dabs(ares(2*i-1)) / xstep
    endif
  enddo
  ndmax = (nmax-ider)/2 + 1

  !------------------------------------------------------------------------------!
  deallocate( expan  )
  deallocate( expan1 )
  deallocate( expans )
  deallocate( ares   )
  deallocate( ycoor2 )

  !------------------------------------------------------------------------------!
  end subroutine centder

  !==============================================================================!
  !==============================================================================!

  subroutine stir1(n,k,istirc1)
  !------------------------------------------------------------------------------!
  !   Stirling numbers of the first kind: [n,k] = (n-1)*[n-1,k] + [n-1,k-1]
  !------------------------------------------------------------------------------!
  implicit none

  integer :: n, k

  integer, dimension(0:n,0:k) :: istirc1

  integer :: ik, jn
  !------------------------------------------------------------------------------!
  istirc1 = 0  ! istirc1(0:nstir)

  do ik=0, k
    istirc1(ik,ik) = 1
  enddo

  do ik=1, k
    do jn=ik+1, n
      istirc1(jn,ik) = (jn-1) * istirc1(jn-1,ik) + istirc1(jn-1,ik-1)
    enddo
  enddo

  !------------------------------------------------------------------------------!
  end subroutine stir1

  !==============================================================================!
  !==============================================================================!

  subroutine mul2exp(nord,n1,n2,a1,a2,ares)
  !------------------------------------------------------------------------------!
  !              Multiply two expansions up tu desired order                     !
  !   It is assumed that both expansions are sorted, so if you want use          !
  !   expansion which have only first, third and sixth term you must send array  !
  !                           (a1, 0, a3, 0, 0, a6)                              !
  !------------------------------------------------------------------------------!
  implicit none

  integer :: nord, n1, n2

  real(8), dimension(0: n1 ) :: a1
  real(8), dimension(0: n2 ) :: a2
  real(8), dimension(0:nord) :: ares

  integer :: i, j
  !------------------------------------------------------------------------------!
  ares = 0.0d0  ! ares(nord)
  do i=0, n1
    if((i+1) == nord)exit
    do j=0, n2
      if((i+j).gt.nord)exit
      ares(i+j) = ares(i+j) + a1(i) * a2(j)
    enddo
  enddo

  !------------------------------------------------------------------------------!
  end subroutine mul2exp

  !==============================================================================!
  !==============================================================================!

  recursive function forwn(N,ycoor,ipoint,npoints,ier) result (resfor)
  implicit none
  ! Calculates recursively N forward differences in point "ipoint"

  integer :: npoints, ipoint, N, ier  ! external variables

  integer :: ip ! internal variables

  real(8) :: resfor
  real(8), dimension(npoints) :: ycoor

  !------------------------------------------------------------------------------!
  ip = ipoint + 1

  if (ip == (npoints+1)) then
    write(luout,'(a)') ' Not enough points to calculate your forward difference !!'
    write(luout,'(a,i4)') ' WARNING in forwn (module num_der); N =',N
    ier = 1; resfor = 0.0d0
  endif

  if ( N==1 ) then 
    resfor = ycoor(ip) - ycoor(ipoint)
  else
    resfor = forwn(N-1,ycoor,ip,npoints,ier) - forwn(N-1,ycoor,ipoint,npoints,ier)
  endif

  !------------------------------------------------------------------------------!
  end function forwn

  !==============================================================================!
  !==============================================================================!

  recursive function backn(N,ycoor,ipoint,npoints,ier) result (resbac)
  implicit none
  ! Calculates recursively N backward differences in point "ipoint"

  integer :: npoints, ipoint, N, ier  ! external variables

  integer :: ip ! internal variables

  real(8) :: resbac
  real(8), dimension(npoints) :: ycoor

  !------------------------------------------------------------------------------!
  ip = ipoint - 1

  if (ip < 0) then
    write(luout,'(a)') ' Not enough points to calculate your backward difference !!'
    write(luout,'(a,i4)') ' WARNING in backn (module num_der); N =',N
    ier = 1; resbac = 0.0d0
  endif

  if ( N==1 ) then 
    resbac = ycoor(ipoint) - ycoor(ip)
  else
    resbac = backn(N-1,ycoor,ipoint,npoints,ier) - backn(N-1,ycoor,ip,npoints,ier)
  endif

  !------------------------------------------------------------------------------!
  end function backn

  !==============================================================================!
  !==============================================================================!

  recursive function centn(N,ycoor,ipoint,npoints,ier) result (resdel)
  implicit none
  ! Calculates recursively N central differences in point "ipoint"

  integer :: npoints, ipoint, N, ier  ! external variables

  integer :: Nm, ip, im ! internal variables

  real(8) :: resdel
  real(8), dimension(npoints) :: ycoor

  !------------------------------------------------------------------------------!
  im = ipoint - 1
  ip = ipoint + 1

  if ( ( im == 0 ) .or. (ip == (npoints+1)) ) then
    write(luout,'(a)') ' Not enough points to calculate your central difference !!'
    write(luout,'(a,i4)') ' WARNING in centn (module num_der); N =',N
    ier = 1; resdel = 0.0d0
  endif

  Nm = N - 1
  if ( N==1 ) then 
    resdel = ycoor(ip) - ycoor(im)
  else
    resdel = centn(Nm,ycoor,ip,npoints,ier) - centn(Nm,ycoor,im,npoints,ier)
  endif

  !------------------------------------------------------------------------------!
  end function centn

  !==============================================================================!
  !==============================================================================!

  recursive function centnmu(N,ycoor,ipoint,npoints,ier) result (resdelmu)
  implicit none
  ! Calculates recursively N central differences in point "ipoint"
  !   and multiply it with averaging operator from the left

  ! external variables
  integer :: npoints, ipoint, N, ier
  real(8), dimension(npoints) :: ycoor

  ! internal variables
  integer :: ip, im
  real(8) :: resdelmu

  !------------------------------------------------------------------------------!
  im = ipoint - 1
  ip = ipoint + 1

  if ( ( im == 0 ) .or. (ip == (npoints+1)) ) then
    write(luout,'(a)') ' Not enough points to calculate your central difference !!'
    write(luout,'(a,i4)') ' WARNING in centnmu (module num_der); N =',N
    ier = 1; resdelmu = 0.0d0
  endif

  if ( N == 0 ) then 
    resdelmu = (ycoor(ip) + ycoor(im)) / 2.0d0
  else
    resdelmu = centnmu(N-1,ycoor,ip,npoints,ier) - centnmu(N-1,ycoor,im,npoints,ier)
  endif

  !------------------------------------------------------------------------------!
  end function centnmu

!================================================================================!
!                             END OF THE MODULE                                  !
!================================================================================!
end module num_der
