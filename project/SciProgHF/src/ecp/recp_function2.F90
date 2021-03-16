!      Copyright (c) 2019 by the authors of DIRAC.
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

MODULE RECP_FUNCTION2
CONTAINS
SUBROUTINE bummer( text, ierr, errtyp )
!  process a program error.
!  input:
!  text  = character string to be printed.
!  ierr  = internal program error to be printed.
!  errtyp = 0 for warning.  traceback may be generated. execution continues.
!         = 1 for nonfatal error.  traceback may be generated.
!             execution is stopped. jcl condition code is set to allow
!             subsequent program steps to continue if possible.
!         = 2 for fatal error.  traceback may be generated.
!             execution is stopped. jcl condition code is set to abort
!             subsequent program steps if possible.
!  entry ibummr must be called prior to bummer() to set the output
!  unit and to perform any additional initialization.
      implicit integer(a-z)
      character*(*) text
      integer iunit, errtyp, ierr
      integer f77err
 
      if ( errtyp .eq. 0 ) then
!        # print a warning message and continue execution.
         write(6,6010) 'bummer (warning):', text, ierr
         return
      elseif ( errtyp .eq. 1 ) then
!        # print a warning message, stop execution.
         write(6,6010) 'bummer (nonfatal):', text, ierr
         call exit( 0 )
      elseif ( errtyp .eq. 2 ) then
!        # print an error message, stop execution, and abort job sequence.
         write(6,6010) 'bummer (fatal):', text, ierr
         call exit( 1 )
      else
!        # unknown error level.  treat as a fatal error.
         write(6,6020) 'bummer (unknown): errtyp=', errtyp, text, ierr
         call exit( 1 )
      endif
 
!     # this statement is not executed, it is included just to avoid compiler warnings. -rls
      stop 'bummer error'
!     # initialization...
      entry ibummr( iunit )
!     # save the listing unit for use later.
      nlist = iunit
      return
6010  format(1x,a,a,i10)
6020  format(1x,a,i10,a,i10)
END SUBROUTINE BUMMER


subroutine siftyp( itypea, itypeb, chrtyp )
!  return a character description of the integral type.
!  input:
!  itypea, itypeb = generic and specific integral or energy(*) types.
!  output:
!  chrtyp = character description. (this should be at least character*8
!           in the calling program.)
      implicit none
      integer          itypea,  itypeb
      character*(*)    chrtyp
!     local variables
      integer  i
      integer  typ1e,  st1e,   end1e
      integer  typ2e,  st2e,   end2e
      integer  typc,   stc,    endc
      integer  typte,  stte,   endte
      integer  typcv,  stcv,   endcv
!     # typ1e = number of defined 1-e array types.
!     # typ2e = number of defined 2-e array types.
!     # typc  = number of core energy types.
!     # typte = number of total energy types.
!     # typcv = number of convergence types.
      parameter( typ1e = 16, st1e = 1,        end1e = typ1e        )
      parameter( typ2e =  2, st2e = end1e +1, end2e = end1e +typ2e )
      parameter( typc  =  1, stc  = end2e +1, endc  = end2e +typc  )
      parameter( typte = 14, stte = endc  +1, endte = endc  +typte )
      parameter( typcv = 10, stcv = endte +1, endcv = endte +typcv )
      integer    ntype
      parameter( ntype = endcv )
      integer ltypea(ntype), ltypeb(ntype)
      character*8 lctype(ntype)
!     # these array types can be in any order, so new ones can be added
!     # to the end or inserted into the middle as appropriate. -rls
!
!     # warning: case dependent code:  do not change the case of the
!     # following character strings.
  data ( ltypea(i), ltypeb(i), lctype(i), i = st1e, end1e ) / &
      0,0,'S1(*)',   0,1,'T1(*)',   0,2,'V1(*)',   0,3,'Veff(*)',  0,4,'VFC(*)', &
      0,5,'Vref(*)', 0,6,'H1(*)',   0,7,'D1(*)',   0,8,'F(*)',     0,9,'Q(*)', &
      1,0,'X(*)',    1,1,'Y(*)',    1,2,'Z(*)',  &
      2,0,'Im(SO:x)',2,1,'Im(SO:y)',2,2,'Im(SO:z)'  / 
      data ( ltypea(i), ltypeb(i), lctype(i), i = st2e, end2e ) / &
      3,0,'1/r12',   3,1,'d2(*)' /
      data ( ltypea(i), ltypeb(i), lctype(i), i = stc, endc ) / &
      0,-1,'Nuc.Rep.' /
      data ( ltypea(i), ltypeb(i), lctype(i), i = stte, endte ) / &
      -1,  0,'SCF',     -1, -1,'MCSCF',  -1, -2,'MRSDCI',   -1, -3,'CPF',     &
      -1, -4,'ACPF',    -1, -5,'LCC-SD', -1, -6,'MRPT',     -1, -7,'Bk',      &
      -1, -8,'DV1',     -1, -9,'DV2',    -1,-10,'EPOPLE',   -1,-11,'S.O. CI', &
      -1,-12,'SR-SDCI', -1,-13,'UCEPA'  /
      data ( ltypea(i), ltypeb(i), lctype(i), i = stcv, endcv ) / &
      -2, 0,'SCF-D.E.', -2,-1,'SCF-D.D1', -2,-2,'MC-D.E.',  -2,-3,'MC-Wnorm', &
      -2,-4,'MC-Knorm', -2,-5,'MC-ApxDE', -2,-6,'Bk-Resid', -2,-7,'CI-Resid', &
      -2,-8,'CI-D.E.',  -2,-9,'CI-ApxDE' /
 
  do i = 1, ntype
     if ( itypea .eq. ltypea(i) ) then
        if ( itypeb .eq. ltypeb(i) ) then
           chrtyp = lctype(i)
           return
        endif
     endif
  enddo
! # loop exit means unrecognized type.
  chrtyp = 'Unknown'
  return
end subroutine siftyp


!----------------------- end colib --------------------

subroutine gcentr(ica,nc)
! generate group of center-interchange operators
  integer        mccu,mconu,mcu,mpru,nbft,nnbft,mrcru,mstu,msu,ng,ns,nst
  common /parmi/ mccu,mconu,mcu,mpru,nbft,nnbft,mrcru,mstu,msu,ng,ns,nst
  dimension ica(mcu,msu,*), nc(*)

      if (ng .eq. 1) goto 56
 
!     set up identity operator
      do 16 is=1,ns
        icu=nc(is)
        if(icu.eq.1) go to 16
        do 12 ic=1,icu
          ica(ic,is,1)=ic
   12   continue
   16 continue
 
!        do 52(+3) ig=2,ng
 
      ig=2
   20 do 52 jg=2,ig
      igp=ig
      jgp=jg
 
!  generate new operator candidate
   24 do 32 is=1,ns
        icu=nc(is)
        if(icu.eq.1) go to 32
        do 28 ic=1,icu
          ica(ic,is,ng+1)=ica(ica(ic,is,igp),is,jgp)
   28   continue
   32 continue
 
!  check against list of operators
      do 44 kg=1,ng
        do 40 is=1,ns
          icu=nc(is)
          if(icu.eq.1) go to 40
          do 36 ic=1,icu
            if(ica(ic,is,ng+1).ne.ica(ic,is,kg)) go to 44
   36     continue
   40   continue
!       # not a new operator
        go to 48
   44 continue
!     # new operator found
      ng=ng+1
   48 if(jgp.eq.ig) go to 52
 
!  take product in other order
      igp=jg
      jgp=ig
      go to 24
   52 continue
      if(ig.eq.ng) go to 56
      ig=ig+1
      go to 20
 
!  set up operators for unique nuclei
   56 do 64 is=1,ns
      if(nc(is).eq.1) then
        do 60 ig=1,ng
          ica(1,is,ig)=1
   60   continue
      endif
   64 continue
      return
end subroutine gcentr


function qcomp(alpha,dfac,n,l,t,xk)
! compute q(n,l) scaled by sqrt(pi)*exp(-t) to prevent overflows
! arguments are alpha, xk, and t=xk**2/(4*alpha)
! no restriction on the magnitude of t
! increase dfac array to raise n, l restrictions
  implicit real*8 (a-h,o-z)
  parameter (am1=-1.0d0, a0=0.0d0, accpow=1.0d-14, &
    accasy=1.0d-10, a1rtpi=0.56418958354775629d0, a1=1.0d0,a2=2.0d0, a4=4.0d0)
  dimension dfac(*)
  dimension tmin(0:8)
  data tmin/31.0d0,28.0d0,25.0d0,23.0d0,22.0d0,20.0d0,19.0d0,18.0d0,15.0d0/
 
  if(mod(n+l,2).ne.0.or.n.le.l) go to 30
 
!  use alternating series (n+l.le.22.and.l.le.10)
      if(l.eq.0) then
        xkp=a1
      else
        xkp=(xk/(alpha+alpha))**l
      endif
      prefac=xkp * dfac(n+l+1)/((alpha+alpha)**((n-l)/2) * sqrt(a4*alpha) * dfac(l+l+3))
      num=l-n+2
      xden=(l+l+3)
      term=a1
      sum=term
      xc=am1
   10 if(num.ne.0) then
        fnum=num
        term=term*fnum*t/(xden*xc)
        xc=xc+am1
        sum=sum+term
        num=num+2
        xden=xden+a2
        go to 10
      endif
 
      qcomp=prefac*sum
      return
 
   30 if(t.lt.tmin(min(n,8))) go to 60
 
!  use asymptotic series (arbitrary n,l)
      xkp=(xk/(alpha+alpha))**(n-2)
      prefac=xkp/((alpha+alpha)*sqrt(a4*alpha))
      sum=a1
      term=a1
      fac1=(l-n+2)
      fac2=(1-l-n)
      xc=a1
   40 term=term*fac1*fac2/(a4*xc*t)
      if(term.eq.a0) go to 50
      sum=sum+term
      if(abs(term/sum).lt.accasy) go to 50
      fac1=fac1+a2
      fac2=fac2+a2
      xc=xc+a1
      go to 40
   50 qcomp=prefac*sum
      return
 
!  use power series (n+l.le.22.and.l.le.10)
   60 if(l.eq.0) then
        xkp=a1
      else
        xkp=(xk/(alpha+alpha))**l
      endif
      prefac=exp(-t)*xkp/(alpha+alpha)**((n-l+1)/2)
      if(mod(n+l,2).eq.0) then
        prefac=prefac/sqrt(a4*alpha)
      else
        prefac=a1rtpi*prefac
      endif
      xnum=(l+n-1)
      xden=(l+l+1)
      term=dfac(l+n+1)/dfac(l+l+3)
      sum=term
      xj=a0
   70 xnum=xnum+a2
      xden=xden+a2
      xj=xj+a1
      term=term*t*xnum/(xj*xden)
      sum=sum+term
      if((term/sum).gt.accpow) go to 70
      qcomp=prefac*sum
      return
end function qcomp

END MODULE RECP_FUNCTION2
