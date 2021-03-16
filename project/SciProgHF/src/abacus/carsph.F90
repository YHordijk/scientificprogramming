!==============================================================================!
module carsph
!------------------------------------------------------------------------------!
!
! NOTE:
!
! (1) Dalton type ordering of spherical basis functions -l,..., 0,...,+l
!
! (2) Contains open-ended transformation matrices 
!
! (3) internal counting of l starts at 1 (not at 0)
!
! (4) normalized to sqrt(l!!) as in Dalton
!
!
! BHP, fall 2016
!
!------------------------------------------------------------------------------!

 implicit none

 private

 logical :: initialized = .false.

 integer :: lmax   = 0
 integer :: ntsph  = 0
 integer :: nttrf  = 0
 integer :: minteg = 0

 integer, pointer :: iofsph(:) => null()
 integer, pointer :: isph(:)   => null()
 integer, pointer :: jcar(:)   => null()

 real(8), pointer :: csph(:)   => null()
 real(8), pointer :: cinv(:)   => null()

 real(8), pointer :: csnrm(:)   => null()
 real(8), pointer :: cinrm(:)   => null()

 public  :: carsph_initialized, carsph_init, carsph_clear, car2sph, sph2car, getnorm

 contains

!==============================================================================!
logical function carsph_initialized()

 implicit none

 carsph_initialized = initialized

end function carsph_initialized
!==============================================================================!
!==============================================================================!
subroutine carsph_init(nhkt,mintg0)
!------------------------------------------------------------------------------!
!
! generate the 
!
!  A) Cartesian -> spherical Harmonic  and
!
!  B) inverse  Cartesian -> spherical Harmonic 
!
! transformation matrices up to maximum angular momentum
!
! compute and store their norm for each angular momentum
!
!------------------------------------------------------------------------------!

 implicit none

! parameter:
 logical, parameter :: locdbg = .false., debug_ortho = .false.
 character(len=*), parameter :: chrdbg = "carsph_init>"

! input:
 integer, intent(in) :: nhkt
 integer, intent(in), optional :: mintg0

! local: 
 integer :: lval, mval, it, iu, iv, nx, ny, nz, nv, nt, ivm, mabs
 integer :: ifct, irow, icol, ideg, irun, ltmp, ldx
 integer :: nsphr, ncart, mxcar, itrf, novl, ioff, ldim, ksph
 integer :: idx, jdx, kdx, lxi, lyi, lzi, lxj, lyj, lzj, ip, jp, kp, nu, nup
 real(8) :: dsgn, dfct, coeff, bilt, bitu, bimv, dnum, deno
 real(8) :: dnorm, dtmp, ttmp, tmin

! local allocatables:
 integer, allocatable :: icart(:), iovl(:), jovl(:)
 integer(8), allocatable :: factorial(:),  doubfact(:)
 real(8), allocatable :: ccart(:), smat(:)
 real(8), allocatable :: scr1(:,:), scr2(:,:), scr3(:,:)

 if (present(mintg0)) then
  minteg = mintg0
 else
  minteg = 2
 end if

 if (minteg.lt.1 .or. minteg.gt.3 ) then
  call quit("minteg should be 1, 2 or 3!")
 end if

 ! flag for initialization
 initialized = .true.

 ! assign maximum angular momentum
 lmax = nhkt-1

 ! determine dimensions
 mxcar = (lmax+1)*(lmax+2)/2

 ntsph = 0
 do lval = 2,lmax
  ntsph = ntsph+2*lval+1
 end do

 ! module pointers
 allocate(iofsph(lmax),isph(ntsph+1),csnrm(lmax-1),cinrm(lmax-1))

 ! temporary arrays
 allocate(icart(mxcar),ccart(mxcar),iovl(mxcar+1),jovl(mxcar*mxcar),&
          smat(mxcar*mxcar),factorial(0:2*lmax),doubfact(-1:2*lmax-1))

 !------------------------------------------------!
 ! 1) determine offset array for the 
 !    spherical Harmonic indeces
 !------------------------------------------------!
 iofsph(1) = 0
 do lval = 2,lmax
  iofsph(lval) = iofsph(lval-1)+2*lval+1
 end do

 itrf = 0
 isph(1) = 1

 !------------------------------------------------!
 ! compute factorial
 !------------------------------------------------!
 factorial(0) = 1
 do lval = 1,2*lmax
  factorial(lval) = factorial(lval-1)*int(lval,8)
 end do

 !------------------------------------------------!
 ! double factorial
 !------------------------------------------------!
 doubfact(-1) = 1
 doubfact( 0) = 1
 doubfact( 1) = 1

 ldim = 2*lmax-1
 ltmp = ldim/2
 ltmp = 2*ltmp
 do lval = 2,ltmp,2
  doubfact(lval) = doubfact(lval-2)*int(lval,8)
 end do

 ltmp= ltmp+mod(ldim,2)
 do lval = 3,ltmp,2
  doubfact(lval) = doubfact(lval-2)*int(lval,8)
 end do

 !------------------------------------------------!
 ! 2a) allocate all arrays
 !  b) compute the Cart -> sph. Harm trafo matrices
 !------------------------------------------------!
 run_loop: do irun = 1,2

  ! allocate arrays for row indeces and coefficients
  if (irun.eq.2) then
   allocate(jcar(nttrf),csph(nttrf),cinv(nttrf))
  end if

  nttrf = 0

  traf_loop: do lval = 2,lmax

   ncart = (lval+1)*(lval+2)/2

   do mval = -lval,lval

    ! initialize with 0
    do ideg = 1,ncart
     icart(ideg) = 0
    end do
    do ideg = 1,ncart
     ccart(ideg) = 0.d0
    end do

    mabs = abs(mval)
    nt = (lval-mabs)/2
    if (mval.ge.0) then
     ivm = 0
     nv = mabs/2
     nv = 2*nv
    else
     ivm = 1
     nv = (mabs-ivm)/2
     nv = 2*nv+ivm
    end if
 
    if (mval.eq.0) then
     deno = sqrt(2.d0)
    else
     deno = 1.d0
    end if
 
    dnum = sqrt(real(2*factorial(lval+mabs)*factorial(lval-mabs),8))
 
    deno = real(factorial(lval)*(2**mabs),8)*deno
 
    ! --------------------------------------- !
    ! extra normalization for minteg = 2 case
    ! --------------------------------------- !
    if (minteg.eq.2 .or. minteg.eq.3) deno = deno * sqrt(real(doubfact(2*lval-1),8))

    dnorm = dnum / deno
 
    do it = 0,nt
     ifct = 4**it
     dfct = 1.d0 / real(ifct,8)
     bilt = binom(factorial,lval,it)*binom(factorial,lval-it,mabs+it)
     do iu = 0,it
      bitu = binom(factorial,it,iu)
      do iv = ivm,nv,2
       bimv = binom(factorial,mabs,iv)
       if (mod(it+(iv-ivm)/2,2).eq.0) then
        dsgn = +1.d0
       else
        dsgn = -1.d0
       end if
 
       coeff = dsgn*dfct*bilt*bitu*bimv
       
       nx = 2*it+mabs-2*iu-iv
       ny = 2*iu+iv
       nz = lval-2*it-mabs
 
       icol = (lval-nx+1)
       irow = nz+1
 
       ideg = irow+icol*(icol-1)/2
       
       icart(ideg) = 1
       ccart(ideg) = ccart(ideg)+dnorm*coeff

       if (locdbg) then
        write(6,"(a,2(i2,1x),ES25.18,3(1x,i2),1x,i3)") chrdbg,lval,mval,dnorm*coeff,nx,ny,nz,ideg
       end if
 
      end do
     end do
    end do
 
    if (irun.eq.1) then
     isph(itrf+1) = nttrf+1
     itrf = itrf+1
    end if

    ! --------------------------------------- !
    ! extra normalization for minteg = 1 case
    ! --------------------------------------- !
    if (minteg.eq.1 .and. irun.eq.2) then
     ! initialize with maximum absolute value
     tmin = 0.d0
     do ideg = 1,ncart
      ttmp = abs(ccart(ideg))
      if (icart(ideg).gt.0 .and. ttmp.gt.1.D-13 .and. ttmp .gt.tmin) then
       tmin = ttmp
       idx = ideg
      end if
     end do
     tmin = ccart(idx)

     ! look for minimum value
     do ideg = 1,ncart
      ttmp = abs(ccart(ideg))
      if (icart(ideg).gt.0 .and. ttmp.gt.1.D-13 .and. ttmp .lt.tmin) tmin = ttmp
     end do
     tmin = 1.d0 / tmin
     do ideg = 1,ncart
      ccart(ideg) = ccart(ideg) * tmin
     end do
    end if 

    if (irun.eq.1) then
     do ideg = 1,ncart
      if (icart(ideg).gt.0 .and. abs(ccart(ideg)).gt.1.D-13) then
       nttrf = nttrf+1
      end if
     end do
    else if (irun.eq.2) then
     do ideg = 1,ncart
      if (icart(ideg).gt.0 .and. abs(ccart(ideg)).gt.1.D-13) then
       nttrf = nttrf+1
       jcar(nttrf) = ideg
       csph(nttrf) = ccart(ideg)
      end if
     end do
    end if

   end do
 
  end do traf_loop

  if (irun.eq.1) then
   isph(itrf+1) = nttrf+1
  end if

 end do run_loop

 ! --------------------- !
 !
 ! flip phase for ... 
 !
 !  l = 3, m_l = -3     (position 7)
 !
 !  l = 4, m_l = +2, -3 (position 5, 7)
 !
 !  in Turbomole minteg = 3 case
 ! --------------------- !
 if (minteg.eq.3) then

  if (lmax.ge.3) then
   !  l = 3, m_l = -3
   do jp = isph(iofsph(2)+1),isph(iofsph(2)+2)-1
    csph(jp) = -csph(jp)
   end do
  end if

  if (lmax.ge.4) then
   !  l = 4, m_l = -3
   do jp = isph(iofsph(3)+2),isph(iofsph(3)+3)-1
    csph(jp) = -csph(jp)
   end do
 
   !  l = 4, m_l = +2
   do jp = isph(iofsph(3)+7),isph(iofsph(3)+8)-1
    csph(jp) = -csph(jp)
   end do
  end if

 end if

 ! ---------------------------------- !
 ! compute the norm of trafo matrices
 ! ---------------------------------- !
 do lval = 3,nhkt
  ldx  = lval-2
  ioff = iofsph(ldx)
  dtmp = 0.d0
  do ksph = 1,2*lval-1
   do jp = isph(ioff+ksph),isph(ioff+ksph+1)-1
    dtmp = dtmp+csph(jp)*csph(jp)
   end do
  end do
  csnrm(ldx) = sqrt(dtmp)
 end do

 if (locdbg) then
  do lval = 3,nhkt
   ldx  = lval-2
!
   write(6,"(/a,i3,1x,ES25.18)") chrdbg//" Spher. -> Cart. trafo matrix:",lval,csnrm(ldx)
!
   ioff = iofsph(ldx)
   do idx = 1,2*lval-1
    do jp = isph(ioff+idx),isph(ioff+idx+1)-1
     jdx = jcar(jp)
     write(6,"(i3,1x,i3,1x,ES25.18)") idx,jdx,csph(jp)
    end do
   end do
  end do
 end if

 !------------------------------------------------!
 ! 3) compute the inverse tranformation matrices
 !------------------------------------------------!

 inv_loop: do lval = 2,lmax

  ncart = (lval+1)*(lval+2)/2
  nsphr = 2*lval+1

  !------------------------------------------------!
  ! 3a) compute Cartesian overlap
  !------------------------------------------------!
  ip = 1

  idx = 0
  do lxi = lval,0,-1
   do lzi = 0,lval-lxi 
    lyi = lval-lxi-lzi
    idx = idx+1

    iovl(idx) = ip

    jdx = 0
    do lxj = lval,0,-1
     do lzj = 0,lval-lxj
      lyj = lval-lxj-lzj
      jdx = jdx+1

      if ((mod(lxi+lxj,2)+mod(lyi+lyj,2)+mod(lzi+lzj,2)) .eq.0) then
       smat(ip) = doubfact(lxi+lxj-1)*doubfact(lyi+lyj-1)*doubfact(lzi+lzj-1)
       jovl(ip) = jdx
       ip = ip+1
      end if 
     end do
    end do

   end do
  end do

  iovl(ncart+1) = ip
  novl = ip-1

  if (locdbg) then
!
   write(6,"(/a,i3)") chrdbg//" Cartesian overlap",lval
!
   do idx = 1,ncart
    do kp = iovl(idx),iovl(idx+1)-1
     jdx = jovl(kp)
     write(6,"(i3,1x,i3,1x,ES25.18)") idx,jdx,smat(kp)
    end do
   end do
  end if

  !------------------------------------------------!
  ! 3b) comute (V^(-1))^T = S * V
  !------------------------------------------------!
  do nu = 1,ncart
   icart(nu) = 0
  end do
 
  ioff = iofsph(lval-1)

  do idx = 1,nsphr
   do jp = isph(ioff+idx),isph(ioff+idx+1)-1
    jdx = jcar(jp)
    do kp = iovl(jdx),iovl(jdx+1)-1
     kdx = jovl(kp)
     if (icart(kdx) .ne. idx) then
      icart(kdx) = idx
      ccart(kdx) = csph(jp) * smat(kp)
     else
      ccart(kdx) = ccart(kdx) + csph(jp) * smat(kp)
     end if
 
    end do
 
   end do
 
   do nup = isph(ioff+idx),isph(ioff+idx+1)-1
    nu = jcar(nup)
    cinv(nup) = ccart(nu)
   end do
 
  end do

 end do inv_loop

 ! ---------------------------------- !
 ! compute the norm of inverse trafo matrices
 ! ---------------------------------- !
 do lval = 3,nhkt
  ldx  = lval-2
  ioff = iofsph(ldx)
  dtmp = 0.d0
  do ksph = 1,2*lval-1
   do jp = isph(ioff+ksph),isph(ioff+ksph+1)-1
    dtmp = dtmp+cinv(jp)*cinv(jp)
   end do
  end do
  cinrm(ldx) = sqrt(dtmp)
 end do

 if (locdbg) then
  do lval = 3,nhkt
   ldx  = lval-2
!
   write(6,"(/a,i3,1x,ES25.18)") chrdbg//" inverse Spher. -> Cart. trafo matrix:",lval,cinrm(ldx)
!
   ioff = iofsph(ldx)
   do idx = 1,2*lval-1
    do jp = isph(ioff+idx),isph(ioff+idx+1)-1
     jdx = jcar(jp)
     write(6,"(i3,1x,i3,1x,ES25.18)") idx,jdx,cinv(jp)
    end do
   end do
  end do
 end if

 if (debug_ortho) then

  do lval = 3,nhkt

!
   write(6,"(/a,i3)") chrdbg//" test:",lval
!
   ldx  = lval-2

   ncart = (lval)*(lval+1)/2
   nsphr = 2*lval-1

   allocate(scr1(nsphr,ncart),scr2(nsphr,ncart),scr3(nsphr,nsphr))

   scr1(:,:) = 0.d0
   scr2(:,:) = 0.d0

   ioff = iofsph(ldx)
   do idx = 1,2*lval-1
    do jp = isph(ioff+idx),isph(ioff+idx+1)-1
     jdx = jcar(jp)
     scr1(idx,jdx) = csph(jp)
     scr2(idx,jdx) = cinv(jp)
    end do
   end do

   call dgemm("n","t",nsphr,nsphr,ncart,&
              1.d0,scr2,nsphr,scr1,nsphr,0.d0,scr3,nsphr)

   call output(scr3,1,nsphr,1,nsphr,nsphr,nsphr,1,6)

   deallocate(scr1,scr2,scr3)

  end do

 end if

 ! release temporary memory
 deallocate(icart,ccart,iovl,jovl,smat,factorial,doubfact)

end subroutine carsph_init
!==============================================================================!
!==============================================================================!
subroutine carsph_clear()
!------------------------------------------------------------------------------!
!
! deallocate and nullify all internal arrays
!
!------------------------------------------------------------------------------!

 implicit none

 initialized = .false.

 lmax  = 0
 ntsph = 0
 nttrf = 0

 if (associated(iofsph)) deallocate(iofsph)
 if (associated(isph))   deallocate(isph)
 if (associated(jcar))   deallocate(jcar)
 if (associated(csph))   deallocate(csph)
 if (associated(cinv))   deallocate(cinv)
 if (associated(csnrm))  deallocate(csnrm)
 if (associated(cinrm))  deallocate(cinrm)
 
 nullify(iofsph,isph,jcar,csph,cinv,csnrm,cinrm) 

end subroutine carsph_clear
!==============================================================================!
real(8) function binom(factorial,nval,kval)

 implicit none

! input:
 integer, intent(in) :: nval, kval
 integer(8), intent(in) :: factorial(0:nval)

 binom = real(factorial(nval),8)/real(factorial(kval)*factorial(nval-kval),8)

end function binom
!==============================================================================!
!==============================================================================!
 subroutine car2sph(sphbas,carbas,iopt,nhkt)
!------------------------------------------------------------------------------!
!
! transform one basis function from Cartesians to spherical Harmonics 
!
!------------------------------------------------------------------------------!

 implicit none

! dimension:
 integer, intent(in) :: nhkt

! input:
 integer, intent(in) :: iopt
 real(8), intent(in) :: carbas(nhkt*(nhkt+1)/2)

! output:
 real(8), intent(out) :: sphbas(2*nhkt-1)

! local scalars:
 character(len=2) :: ctmp
 integer :: ksph, kcar, ioff, jp
 real(8) :: sphk

! local pointer:
 real(8), pointer :: cpnt(:)

 select case (iopt)
  case (1)
   cpnt => csph
  case (2)
   cpnt => cinv
  case default
   call quit("Only iopt = 1,2 valid!")
 end select

 if (nhkt.gt.lmax+1) then
  write(ctmp,"(i2)") nhkt
  call quit("Angular momentum number "//ctmp//" too high!")
 else if (nhkt.eq.1) then
  sphbas(1) = carbas(1)
 else if (nhkt.eq.2) then
  sphbas(1) = carbas(1)
  sphbas(2) = carbas(2)
  sphbas(3) = carbas(3)
 else
  ioff = iofsph(nhkt-2)
  do ksph = 1,2*nhkt-1
   sphk = 0.d0
   do jp = isph(ioff+ksph),isph(ioff+ksph+1)-1
    kcar = jcar(jp)
    sphk = sphk+cpnt(jp)*carbas(kcar)
   end do
   sphbas(ksph) = sphk
  end do
 end if

 end subroutine car2sph
!==============================================================================!
!==============================================================================!
 subroutine sph2car(carbas,sphbas,iopt,nhkt)
!------------------------------------------------------------------------------!
!
! transform one basis function from spherical Harmonics to Cartesians
!
!------------------------------------------------------------------------------!

 implicit none

! dimension:
 integer, intent(in) :: nhkt

! input:
 integer, intent(in) :: iopt
 real(8), intent(in) :: sphbas(2*nhkt-1)

! output:
 real(8), intent(out) :: carbas(nhkt*(nhkt+1)/2)

! local scalars:
 character(len=2) :: ctmp
 integer :: ksph, kcar, ioff, jp
 real(8) :: sphk

! local pointer:
 real(8), pointer :: cpnt(:)

 select case (iopt)
  case (1)
   cpnt => csph
  case (2)
   cpnt => cinv
  case default
   call quit("Only iopt = 1,2 valid!")
 end select

 if (nhkt.gt.lmax+1) then
  write(ctmp,"(i2)") nhkt
  call quit("Angular momentum number "//ctmp//" too high!")
 else if (nhkt.eq.1) then
  carbas(1) = sphbas(1)
 else if (nhkt.eq.2) then
  carbas(1) = sphbas(1)
  carbas(2) = sphbas(2)
  carbas(3) = sphbas(3)
 else
  do kcar = 1,nhkt*(nhkt+1)/2
   carbas(kcar) = 0.d0
  end do
  ioff = iofsph(nhkt-2)
  do ksph = 1,2*nhkt-1
   sphk = sphbas(ksph)
   do jp = isph(ioff+ksph),isph(ioff+ksph+1)-1
    kcar = jcar(jp)
    carbas(kcar) = carbas(kcar)+cpnt(jp)*sphk
   end do
  end do
 end if

 end subroutine sph2car
!==============================================================================!
!==============================================================================!
real(8) function getnorm(nhkt,iopt)

 implicit none

! input:
 integer,intent(in) :: nhkt, iopt

! local scalars:
 character(len=2) :: ctmp

! local pointer:
 real(8), pointer :: cpnt(:)

 select case (iopt)
  case (1)
   cpnt => csnrm
  case (2)
   cpnt => cinrm
  case default
   call quit("Only iopt = 1,2 valid!")
 end select

 if      (nhkt.gt.lmax+1) then
  write(ctmp,"(i2)") nhkt
  call quit("Angular momentum number "//ctmp//" too high!")
 else if (nhkt.eq.1) then
  getnorm = 1.d0
 else if (nhkt.eq.2) then
  getnorm = 1.d0
 else if (nhkt.ge.3) then
  getnorm = cpnt(nhkt-2)
 end if

end function getnorm
!==============================================================================!

end module carsph
!==============================================================================!
