module second_order_minimization

use func_PipekMezey
implicit none

!=============================================================================================!
! func_PipekMezey:                                                                            !
!   module for calculation of Dirac Pipek-Mezey localization functional                       !
!   based either on Mulliken population or projection analysis                                !
!   --> S. Dubillard, J.-B. Rota, and T. Saue, J.Chem.Phys. 124, 154307 (2006)                !
!   --> http://jcp.aip.org/resource/1/jcpsa6/v124/i15/p154307_s1                              !
!                                                                                             !
!   This subroutines (functions) must be provided by some module, in this case                !
!   it is func_PipekMezey:                                                                    !
!   init_func    --> constructor                                                              !
!   quit_func    --> destructor                                                               !
!   mv2point     --> Move to the new expansion point                                          !
!                --> IN: rejected --> the step was rejected, go to the previous point         !
!                --> IN: restore  --> restore previous point                                  !
!                --> IN: xvec     --> new point in the examined space                         !
!   get_func     --> function: Get value of your functional in current point                  !
!   get_gradient --> Get gradient of your functional in current point                         !
!                --> OUT: gradient                                                            !
!                --> OUT: dnrmg --> norm of the gradient                                      !
!   get_hessian  --> Get Hessian of your functional in the current point                      !
!                --> IN : full --> True : calculate full hessian                              !
!                              --> False: calculate diagonal approximation of the Hessian     !
!                --> OUT: hmat --> Hessian or diagonal approximation of the Hessian           !
!                --> OUT: totdiag --> sum of absolute values of diagonal elements             !
!                --> OUT: totnond --> sum of absolute values of non-diagonal elements         !
!   get_date     --> Get date                                                                 !
!=============================================================================================!
! Input variables of the second_order_minimization module:                                    !
!                                                                                             !
!   ndim    - dimension of the examined space                                                 !
!   iprint  - print level (no prints =0)                                                      !
!   maxitr  - maximum number of iterations                                                    !
!   l       - print unit (usually =6)                                                         !
!                                                                                             !
!   ---- Convergence keywords ----                                                            !
!                                                                                             !
!   heskey = 'DIAG' - diagonal approximation of the hessian is used                           !
!          = 'FULL' - calculation of the full hessian is used                                 !
!          = 'COMB' - there are two stages:                                                   !
!               1: in the first stage only diagonal approximation of the hessian is used      !
!               2: in the second stage the calculation of the full hessian is used            !
!   lgfull - it is used when heskey = FULL or in the second stage of the heskey = COMB option !
!            True  - converge when the difference between the functional value                !
!                    of two subsequent iterations is lower than threshold (thfull)            !
!            False - the gradient criterion or criterion of the number                        !
!                    of negative eigenvalues equals zero will be used                         !
!   thfull - threshold; see lgfull for the explanation                                        !
!   lggrad - it is used when heskey = FULL or heskey = DIAG or                                !
!            in the second stage of the heskey = COMB option                                  !
!            True  - converge when the difference between the gradient of the functional      !
!                    of two subsequent iterations is lower than threshold (thgrad)            !
!            False - the functional value criterion or criterion of the number                !
!                    of negative eigenvalues equals zero will be used                         !
!   thgrad - threshold; see lggrad for the explanation                                        !
!   lgdiag - if heskey = DIAG it have the the same purpose as lgfull;                         !
!            for threshold is used thdiag variable                                            !
!          - if heskey = COMB it is used to switch from the first to the second               !
!            stage of the convergence; it will switch                                         !
!            True  - when the difference between the functional value                         !
!                    of two subsequent iterations is lower than threshold (thdiag)            !
!            False - if the number of negative eigenvalues equals zero                        !
!   thdiag - threshold; see lgdiag for the explanation                                        !
!   lgchck - works with heskey = COMB option                                                  !
!            True  - only one calculation of the full hessian will be performed               !
!=============================================================================================!

public  :: minimization_drv

private :: init_minimization, titler, exp_func, get_X, hdiag, get_mu, updateT

type :: minimization_flags
  character*4 :: heskey
  logical     :: lgfull, lgdiag, lggrad, lgchck
  integer     :: iprint, maxitr, l, ndim
  real*8      :: thfull, thdiag, thgrad
end type minimization_flags

type( minimization_flags ), private :: mnf

! global parameters
real*8, parameter, private :: thresh = 1.0d-13, trsh_mu = 1.0d-4
real*8, parameter, private :: d0 = 0.0d0, d1 = 1.0d0, d2 = 2.0d0, d4 = 4.0d0

contains

  subroutine init_minimization
  !==============================================================================!
  !               Initialize minimization for Dirac program                      !
  !==============================================================================!
  implicit none

  ! Internal variables
  integer     :: ier, iunit
  character*8 :: chdum

  !------------------------------------------------------------------------------!
  !-----------------------------!
  ! Read inp data from the file !
  !-----------------------------!
  iunit = 60

  ! open file
  open(iunit,file='OPTIM.INP',status='old',form='formatted',access='sequential',iostat=ier)

  if(ier /= 0)then
    write(6,*) 'init_minimization: ERROR while opening OPTIM.INP file'
    stop
  endif

  rewind 60

  read(iunit,*) chdum,mnf%iprint
  read(iunit,*) chdum,mnf%maxitr
  read(iunit,*) chdum,mnf%l
  read(iunit,*) chdum,mnf%thfull
  read(iunit,*) chdum,mnf%thdiag
  read(iunit,*) chdum,mnf%thgrad
  read(iunit,*) chdum,mnf%lgfull
  read(iunit,*) chdum,mnf%lgdiag
  read(iunit,*) chdum,mnf%lggrad
  read(iunit,*) chdum,mnf%lgchck
  read(iunit,*) chdum,mnf%heskey

  ! close file
  close(iunit,status='keep')

  !-----------------------------------!
  ! initialize module func_PipekMezey !
  !-----------------------------------!
  ! mnf%ndim  -->  dimension of your space
  call init_func(mnf%ndim)

  end subroutine init_minimization

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine minimization_drv(unit_test)
  !==============================================================================!
  !                                                                              !
  !                  Driver for second-order minimization                        !
  !                                                                              !
  !       funcg --> current functional value                                     !
  !       funco --> functional value from previous iteration                     !
  !       funcq --> quadratic approximation to the functional value              !
  !       xvec  --> point (X) in studied vector space                            !
  !       xvect --> transformed point (U^T*X) in studied vector space            !
  !                                                                              !
  !==============================================================================!
  implicit none

  ! External variables
  logical, optional :: unit_test

  ! Internal variables
  character :: daytid*24, hessian*4, chstart*23, resst*50
  logical   :: full, rejected, Iminmin, trdef, nstep, restore

  integer :: itr, ist, ier, nep, nem, nez, nref

  real*8  :: funcg, funcq, funco, dnrmk, dnrmg, diffg, diffq, diffa
  real*8  :: totdiag, totnond, shift, tradius, tpred

  real*8,  allocatable, dimension(:) :: grad, etmp, xvec

  integer, pointer :: iorig(:)
  real*8,  pointer :: xvect(:), gradt(:), eval(:), hmat(:,:)

  !------------------------------------------------------------------------------!
  ! initialaze minimization module
  call init_minimization

  ! Check thresholds
  if(mnf%thdiag < thresh)then
    write(6,'(/a)')      ' WARNING: The threshold (THDIAG) from the input is too small '
    write(6,'(a,d8.1)')  ' WARNING: It will be changed to the minimal reasonable value: ',thresh
    mnf%thdiag = thresh
  endif
  if(mnf%thfull < thresh)then
    write(6,'(/a)')      ' WARNING: The threshold (THFULL) from the input is too small '
    write(6,'(a,d8.1)')  ' WARNING: It will be changed to the minimal reasonable value: ',thresh
    mnf%thfull = thresh
  endif

  ! Set defaults
  full = .false.
  if(mnf%heskey == 'FULL') full = .true.
  hessian = 'DIAG'
  if(full) hessian = 'FULL'

  rejected = .false.
  Iminmin  = .false.
  chstart  = 'START ITERATION NO.    '

  tradius = 10.0d0
  trdef   = .true.
  nstep   = .true.
  restore = .false.

  if(mnf%iprint.ge.4)then
    write(mnf%l,*) 'full  =',full
    write(mnf%l,*) 'thdiag=',mnf%thdiag
    write(mnf%l,*) 'thfull=',mnf%thfull
    write(mnf%l,*) 'thgrad=',mnf%thgrad
    write(mnf%l,*) 'lgdiag=',mnf%lgdiag
    write(mnf%l,*) 'lgfull=',mnf%lgfull
    write(mnf%l,*) 'lggrad=',mnf%lggrad
    write(mnf%l,*) 'lgchck=',mnf%lgchck
    write(mnf%l,*) 'heskey=',mnf%heskey
  endif

  ! Memory allocation
  allocate( grad(mnf%ndim), stat=ier ); call errall(ier,'minimization_drv','grad')
  allocate( etmp(mnf%ndim), stat=ier ); call errall(ier,'minimization_drv','etmp')
  allocate( xvec(mnf%ndim), stat=ier ); call errall(ier,'minimization_drv','xvec')
  grad = d0; etmp = d0; xvec = d0

  ! nullify pointers for later allocation of pointers
  nullify(iorig); nullify(gradt); nullify(eval); nullify(xvect); nullify(hmat)

  !------------------------------------------------------------------------------!
  !                           Main iteration procedure                           !
  !------------------------------------------------------------------------------!
  call titler(' REFERENCE POINT ','=',122)

  funcg = get_func()
  funco = funcg

  write(mnf%l,'(/a,2x,a4, 4x,a, 6x,a, 5x,a, 4x,a, 3x,a, 4x,a, 5x,a, 6x,a)') &
    '*T','ITR','T_old','T_new ','|KAPPA|','  R   ','REJ','CON','DIFFG','DIFFQ'
  write(mnf%l,'(/a,2x,a4, 4x,a, 6x,a, 5x,a, 2x,a, 3x,a, 5x,a, 6x,a, 6x,a)') &
    '*C','ITR','DIFFA','|GRAD|','|KAPPA|','H(nonD/D)','e: +,-,0','H','DIFFG','DIFFQ'

  ! zero step
  call get_hessian(hmat,full,totdiag,totnond)
  call get_gradient(grad,dnrmg)
  call hdiag(hmat,grad,etmp,gradt,eval,iorig,nep,nem,nez,Iminmin,full)

  shift = d0
  dnrmk = d0
  funcq = d0; diffq = d0
  funcg = d0; diffg = d0
  diffa = d0
  itr   = 0
  ist   = 0
  write(mnf%l,'(a,2x,i4,3d11.2,f10.5,1x,3i4,3x,a,2d11.2)') &
    '*C',itr,diffa,dnrmg,dnrmk,totnond/totdiag,nep,nem,nez,hessian,diffg,diffq

  do
    ist = ist + 1

    ! start i-th iteration
    if( .not.rejected .and. .not.restore )then
      write(chstart(20:23),'(i4)') itr+1
      call titler(chstart,'=',122)
      call get_date(daytid)
      write(mnf%l,'(/6x,a,a24)')'--- ',daytid
    endif

    ! get new kappa
    shift = d0
    if(rejected) call get_mu(eval,gradt,tradius,nref,shift,nstep,tpred)
    call get_X(xvect,xvec,gradt,eval,hmat,iorig,nstep,full,tradius,shift,dnrmk)
    funcq = exp_func(xvect,gradt,eval)

    if( (rejected) .and. (dabs(tpred-dnrmk) > trsh_mu*tradius) )then
      write(mnf%l,'(a)') ' WARNING: Predicted and actual norm of kappa are too different'
      write(mnf%l,'(a,i4,a,d12.4,a,d12.4,a,2d12.4)')&
        ' WARNING: itr:',itr,' Predicted:',tpred,' Actual:',dnrmk,' Trust-radius:',tradius
    endif

    ! move to the new expansion point
    call mv2point(restore,rejected,xvec)

    ! get functional value in the new point (new MO)
    funcg = get_func()

    ! Update Trust-radius
    diffg = funcg - funco
    diffq = funcq
    call updateT(dnrmk,diffg,diffq,tradius,rejected,trdef,nstep,restore,full,itr,nem,&
                 gradt,eval,nref,dnrmg)

    ! Move to the new point if the step was not rejected or old step is going to be restored
    if( .not.rejected .and. .not.restore )then
      !---------------------------------------------!
      ! Calculate hessian and gradient in new point !
      !---------------------------------------------!
      call get_hessian(hmat,full,totdiag,totnond)
      call get_gradient(grad,dnrmg)
      call hdiag(hmat,grad,etmp,gradt,eval,iorig,nep,nem,nez,Iminmin,full)

      itr = itr + 1

      if((diffg > d0).and.(diffg > thresh)) &
        write(mnf%l,'(a,i4)')' WARNING: dinvm_new - dinvm_old > 0; itr =',itr
      diffa = dabs(funcg - funco)
      hessian = 'DIAG'
      if(full) hessian = 'FULL'
      funco = funcg

      write(mnf%l,'(/a,2x,a4, 4x,a, 6x,a, 5x,a, 2x,a, 3x,a, 5x,a, 6x,a, 6x,a)') &
        '  ','ITR','DIFFA','|GRAD|','|KAPPA|','H(nonD/D)','e: +,-,0','H','DIFFG','DIFFQ'
      write(mnf%l,'(a,2x,i4,3d11.2,f10.5,1x,3i4,3x,a,2d11.2)') &
        '*C',itr,diffa,dnrmg,dnrmk,totnond/totdiag,nep,nem,nez,hessian,diffg,diffq

      !----------------------------!
      ! check convergence criteria !
      !----------------------------!
      if(mnf%heskey == 'DIAG')then
        if( .not.mnf%lggrad .and.   .not.mnf%lgdiag    .and. Iminmin ) exit
        if(      mnf%lggrad .and. (dnrmg < mnf%thgrad) .and. Iminmin ) exit
        if(      mnf%lgdiag .and. (diffa < mnf%thdiag) .and. Iminmin ) exit
      endif

      if(mnf%heskey == 'FULL')then
        if( .not.mnf%lggrad .and. .not.mnf%lgfull      .and. Iminmin ) exit
        if(      mnf%lggrad .and. (dnrmg < mnf%thgrad) .and. Iminmin ) exit
        if(      mnf%lgfull .and. (diffa < mnf%thfull) .and. Iminmin ) exit
      endif

      if(mnf%heskey == 'COMB')then
        if(full)then
          if( .not.mnf%lggrad .and. .not.mnf%lgfull      .and. Iminmin ) exit
          if(      mnf%lggrad .and. (dnrmg < mnf%thgrad) .and. Iminmin ) exit
          if(      mnf%lgfull .and. (diffa < mnf%thfull) .and. Iminmin ) exit
        else
          if( .not.mnf%lgdiag                            .and. Iminmin ) full = .true.
          if(      mnf%lgdiag .and. (dnrmg < mnf%thdiag) .and. Iminmin ) full = .true.
          if(full)then
            call get_hessian(hmat,full,totdiag,totnond)
            call get_gradient(grad,dnrmg)
            call hdiag(hmat,grad,etmp,gradt,eval,iorig,nep,nem,nez,Iminmin,full)
            hessian = 'FULL'
            write(mnf%l,'(a,2x,i4,3d11.2,f10.5,1x,3i4,3x,a,2d11.2)') &
              '*C',itr,diffa,dnrmg,dnrmk,totnond/totdiag,nep,nem,nez,hessian,diffg,diffq
            if(mnf%lgchck) exit
          endif
        endif
      endif

      if(itr == mnf%maxitr)exit
    endif

    if(ist > 2*mnf%maxitr)exit

  enddo

  if(Iminmin)then
    write(resst,'(a,i4,a,i4,a)') 'MAXIMUM reached after',itr,' iterations and',ist,' steps'
  elseif( full .and. ((diffa < mnf%thfull).and.(dnrmg < mnf%thgrad)) )then
    write(resst,'(a,i4,a,i4,a)') '   Convergence in',itr,' iterations and',ist,' steps    '
  else
    write(resst,'(a,i4,a,i4,a)') 'No convergence after',itr,' iterations and',ist,' steps '
  endif
  call titler(resst,'*',104)

  !------------------------------------------------------------------------------!
  !                Write information for unit test to the disk                   !
  !------------------------------------------------------------------------------!
  if(present(unit_test))then
    ! open file
    open(60,file='OPTIM.OUT',status='unknown',form='formatted',access='sequential',iostat=ier)
    if(ier /= 0)then
      write(6,*) 'minimization_drv: ERROR while creating new OPTIM.OUT file'
      write(0,*) 'minimization_drv: ERROR while creating new OPTIM.OUT file'
      stop 1
    endif
    rewind 60

    write(60,*) funcg

    ! close file
    close(60,status='keep')
  endif

  !------------------------------------------------------------------------------!
  !                        Deallocate arrays or pointers                         !
  !------------------------------------------------------------------------------!
  call quit_func(hmat)

  deallocate( grad )
  deallocate( etmp )
  deallocate( xvec )

  if(associated(iorig)) deallocate( iorig )
  if(associated(xvect)) deallocate( xvect )
  if(associated(gradt)) deallocate( gradt )
  if(associated(eval )) deallocate( eval  )

  end subroutine minimization_drv

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine titler(head,st,in)
  implicit none

  ! External variables
  character, intent(in) :: head*(*), st*(*)
  integer,   intent(in) :: in

  ! Internal variables
  integer   :: i, lhead, length, marg, indent

  !------------------------------------------------------------------------------!
  lhead  = len(head)
  length = lhead

  marg = (80 - length)/2 - 1
  marg = min(marg, in - 100)
  if (marg .gt. 0) marg = marg + 1
  length = length + 2*marg

  indent = max(1,(80 - length)/2 + 1)

  write (mnf%l,'(//150a)') (' ',i=1,indent),(st,i=1,length)
  write (mnf%l,'(150a)') (' ',i=1,indent),(st,i=1,marg-1),' ',head(1:lhead),' ',(st,i=1,marg-1)
  write (mnf%l, '(150a)') (' ',i=1,indent),(st,i=1,length)
  write (mnf%l, '()')

  end subroutine titler

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

  end subroutine

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  real*8 function exp_func(xvect,gradt,eval)
  !==============================================================================!
  !       Calculate value of the quadratic approximation to the functional       !
  !==============================================================================!
  implicit none

  ! External variables
  real*8,  intent(in), pointer :: xvect(:), gradt(:), eval(:)

  ! Internal variables
  integer :: i
  real*8  :: funcq

  !------------------------------------------------------------------------------!
  funcq = d0
  do i=1,size(eval,1)
    funcq = funcq + d2*xvect(i)*gradt(i) + xvect(i)*eval(i)*xvect(i)
  enddo
  exp_func = funcq

  end function exp_func

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine get_X(xvect,xvec,gradt,eval,hmat,iorig,nstep,full,tradius,shift,dnrmk)
  !******************************************************************************!
  !                                                                              !
  !             Construction of the X vector using full hessian                  !
  !          The linear set of equations is solved using information             !
  !                   from diagonalization of the hessian                        !
  !                                                                              !
  !   (H-shift) * X = -g     -->     (tilde{H}-shift) * tilde{X} = - tilde{g}    !
  !                                                                              !
  !   tilde{H} = U^T * H * U = eval                 U = hmat (iorig)             !
  !   tilde{X} = U^T * X     = xvect                X = xvec                     !
  !   tilde{g} = U^T * g     = gradt                dnrmk = |X|                  !
  !                                                                              !
  !        Note: Be careful about minus in imaginary part of gradient in         !
  !                                H * K = - g                                   !
  !                                 H                     K            g         !
  !  / A   B  \   -->   / A^R + B^I   B^I - A^I \      / K^R \      / g^R \      !
  !  \ B*  A* /   -->   \ B^I + A^I   A^R - B^I /      \ K^I /      \-g^I /      !
  !                                                                              !
  !******************************************************************************!
  implicit none

  ! External variables
  logical, intent(in   ) :: full
  logical, intent(inout) :: nstep
  real*8,  intent(in   ) :: shift, tradius
  real*8,  intent(  out) :: dnrmk

  real*8, intent(in ) :: hmat(:,:)
  real*8, intent(out) :: xvec(:)

  integer, intent(in ), pointer :: iorig(:)
  real*8,  intent(in ), pointer :: gradt(:), eval(:)
  real*8,  intent(out), pointer :: xvect(:)

  ! Internal variables
  integer :: i, j, ndim, nred, ier
  real*8  :: dnrm2, dtmp

  !------------------------------------------------------------------------------!
  ! Get some dimensions
  ndim  = size(hmat,1)
  nred  = size(eval,1)

  ! allocate arrays
  if(associated(xvect)) deallocate( xvect )
  allocate( xvect(nred), stat=ier ); call errall(ier,'get_X','xvect(pointer)')

  ! calculate tilde{X}
  if(nstep)then
    do i=1,nred
      xvect(i) = - gradt(i) / (eval(i) - shift)
    enddo
  else
    dtmp = d0
    do i=1,nred
      if(eval(i) > d0) exit
      dtmp = dtmp + eval(i)*eval(i)
    enddo
    dtmp = dsqrt(dtmp)  ! norm of the negative eigenvalues
    ! make step in the direction of negative eigenvalues
    ! eigenvalues are weighted with their relative size 
    xvect(:) = d0
    do i=1,nred
      if(eval(i) > d0) exit
      xvect(i) = tradius * dabs(eval(i)) / dtmp
    enddo
  endif

  ! print \tilde{X}
  if(mnf%iprint.ge.1)then
    write(mnf%l,'(/a)') ' U^T*X vector: '
    write(mnf%l,'(8d12.4)') (xvect(i),i=1,nred)
  endif

  ! calculate X
  xvec(:) = d0
  do j=1,nred
    if(full)then
      do i=1,ndim
        xvec(i) = xvec(i) + hmat(i,j)*xvect(j)
      enddo
    else
      xvec(iorig(j)) = xvect(j)
    endif
  enddo

  ! Norm of the X
  dnrmk = dnrm2(ndim,xvec,1)
  if(mnf%iprint.ge.1)then
    write(mnf%l,'(/a)') ' X vector: '
    write(mnf%l,'(8d12.4)') (xvec(i),i=1,ndim)
    write(mnf%l,'(/a,e12.2)') ' Norm of the kappa: ',dnrmk
  endif

  end subroutine get_X

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine hdiag(hmat,gradient,etmp,gradt,eval,iorig,nep,nem,nez,Iminmin,full)
  !******************************************************************************!
  !                                                                              !
  !               Diagonalize real representation of the hessian                 !
  !                                                                              !
  !******************************************************************************!
  implicit none

  ! External variables
  logical, intent(in ) :: full
  logical, intent(out) :: Iminmin
  integer, intent(out) :: nep, nem, nez
  real*8,  intent(in ), dimension(:) :: gradient
  real*8,  intent(out), dimension(:) :: etmp

  integer, intent(  out), pointer :: iorig(:)
  real*8,  intent(inout), pointer :: hmat(:,:)
  real*8,  intent(  out), pointer :: gradt(:), eval(:)

  ! Internal variables
  integer :: i, j, k, ndim, nred, ier, lwork, ilaenv, nb, icur
  real*8  :: acur

  real*8, allocatable, dimension(:) :: work

  !------------------------------------------------------------------------------!
  ! Get some dimensions
  ndim = size(hmat,1)

  if(full)then

    nb = ilaenv(1,'dsytrd','u',ndim,-1,-1,-1)
    lwork = (nb+2)*ndim  ! see the comment in dsyev lapack routine
  
    allocate( work(lwork), stat=ier ); call errall(ier,'hdiag','work')
  
    ! get eigenvalues (etmp) and eigenvectors (hmat)
    call dsyev('V','U',ndim,hmat,ndim,etmp,work,lwork,ier)
    if(ier /= 0)then
      write(mnf%l,*) ' ERROR: The dsyev LAPACK subroutine in hdiag finished unsuccessfully'
      write(mnf%l,*) ' info =',ier
      stop
    endif

    deallocate( work )

  else
  
    if(associated(iorig)) deallocate( iorig )
    allocate( iorig(ndim), stat=ier ); call errall(ier,'hdiag','iorig(pointer)')

    ! get eigenvalues (etmp)
    ! We use diagonal approximation so it will just sort the diagonal of the hessian
    ! iorig will contain information how the eigenvalues were sorted
    do i=1,ndim
      etmp(i) = hmat(i,1)
      iorig(i)= i
    enddo
    do i=1,ndim
      acur = etmp(i)
      icur = i
      do j=i+1,ndim
        if(acur > etmp(j))then
          acur = etmp(j)
          icur = j
        endif
      enddo
      etmp(icur)  = etmp(i)
      etmp(i)     = acur
      j = iorig(icur)
      iorig(icur) = iorig(i)
      iorig(i)    = j
    enddo

  endif

  ! allocate unnamed arrays (dimension is number of eigenvalues over some threshold)
  nred = 0
  do i=1,ndim
    if(dabs(etmp(i)) > thresh) nred = nred + 1
  enddo
  if(associated(eval))then
    deallocate( gradt )
    deallocate( eval  )
  endif
  allocate( gradt(nred), stat=ier ); call errall(ier,'hdiag','gradt(pointer)')
  allocate( eval (nred), stat=ier ); call errall(ier,'hdiag','eval(pointer)')

  ! project out redundant eigevalues and corresponding eigenvectors
  k = 0
  do i=1,ndim
    if(dabs(etmp(i)) > thresh)then
      k = k + 1
      eval(k) = etmp(i)
      if(full)then
        do j=1,ndim
          hmat(j,k) = hmat(j,i)
        enddo
      else
        iorig(k) = iorig(i)
      endif
    endif
  enddo

  ! calculate transformed gradient: tilde{g} = U^T * g
  gradt(:) = d0
  do i=1,nred
    if(full)then
      do j=1,ndim
        gradt(i) = gradt(i) + hmat(j,i) * gradient(j)
      enddo
    else
      gradt(i) = gradient(iorig(i))
    endif
  enddo

  ! count number of positive, negative and zero eigenvalues
  nep=0; nem=0; nez=0
  do i=1,ndim
    if(dabs(etmp(i)) > thresh)then
      if(etmp(i) < d0)nem = nem + 1
      if(etmp(i) > d0)nep = nep + 1
    else
      nez = nez + 1
    endif
  enddo
  if(nem == 0)Iminmin = .true.  ! we have reached the minimum of the function
  if(nem /= 0)Iminmin = .false.

  if(mnf%iprint.ge.1)then
    write(mnf%l,'(/a,3i4)')' Eigenvalues of the hessian;  +,-,0  =  ',nep,nem,nez
    write(mnf%l,'(8d12.4)') (etmp(i),i=1,ndim)
  endif

  if(mnf%iprint.ge.4)then
    write(mnf%l,'(/a)')' Reduced eigenvalues of the hessian:'
    write(mnf%l,'(8d12.4)') (eval(i),i=1,nred)
    write(mnf%l,'(/a)')' Transformed gradient:'
    write(mnf%l,'(8d12.4)') (gradt(i),i=1,nred)
    write(mnf%l,'(/a)')' Eigenvectors:'
    do i=1,nred
      write(mnf%l,'(a,i5)') ' Vector:',i
      write(mnf%l,'(8d12.4)') (hmat(j,i),j=1,ndim)
    enddo
  endif

  end subroutine hdiag

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine get_mu(eval,gradt,tradius,nref,shift,nstep,tpred)
  !******************************************************************************!
  !                                                                              !
  !                   Calculate level-shift parameter mu                         !
  !                                                                              !
  !          tradius = sqrt{ sum{ (gradt)^2_i / (eval_i-shift)^2  }  }           !
  !                                                                              !
  !******************************************************************************!
  implicit none

  ! External variables
  logical, intent(in ) :: nstep
  integer, intent(in ) :: nref
  real*8,  intent(in ) :: tradius
  real*8,  intent(out) :: shift, tpred

  real*8, intent(in), pointer :: eval(:), gradt(:)

  ! Internal variables
  integer :: i, ndim, imax
  real*8  :: tcurr, xcurr, cutoff, xmax, xmin, step

  !------------------------------------------------------------------------------!
  ! If we are not performing Newton step quit.
  if(.not.nstep)then
    tpred = tradius
    return
  endif

  ! Get some dimensions
  ndim = size(eval,1)

  ! initialization
  imax  = 10000
  cutoff= trsh_mu * tradius
  step  = 0.1d0
  xmax  = eval(nref)
  xmin  = eval(nref) - step
  tcurr = f(xmin)

  ! find lower boundary if necessary
  i = 0
  if(tcurr > tradius)then
    do
      i = i + 1
      if(i > imax)then
        write(mnf%l,'(a)')' WARNING: Maximum number of iteration reached in get_mu(1)!!!'
        write(mnf%l,'(10x,a,d12.4,a,d12.4)')'f(xmin) =',tcurr,' tradius =',tradius
        exit
      endif
      step = xmax - xmin
      xmax = xmin
      xmin = xmin - d2*step
      tcurr = f(xmin)
      if(tcurr < tradius)exit
    enddo
  endif 

  i = 0
  do
  i = i + 1
  if(i > imax)then
    write(mnf%l,'(a)')' WARNING: Maximum number of iteration reached in get_mu(2)!!!'
    write(mnf%l,'(10x,a,d12.4,a,d12.4)')&
      'h(desired) =',tradius,' |h(desired)-h(current)| =',dabs(tradius - tcurr)
    exit
  endif

  if(dabs(tcurr-tradius) < cutoff) exit

  xcurr = (xmax+xmin)/d2
  tcurr = f(xcurr)

  if(tcurr > tradius) xmax = xcurr
  if(tcurr < tradius) xmin = xcurr

  enddo

  shift = xcurr
  tpred = f(xcurr)

  !------------------------------------------------------------------------------!

  contains

    function f(x)
    implicit none
    real*8, intent(in ) :: x
    real*8  :: f
    integer :: ii

    f = d0
    do ii=1,ndim
      if(dabs(gradt(ii)) > thresh)then
        f = f + gradt(ii) * gradt(ii) / ((eval(ii) - x) * (eval(ii) - x))
      endif
    enddo
    f = dsqrt(f)

    end function f

  end subroutine get_mu

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine updateT(dnrmk,diffg,diffq,tradius,rejected,trdef,nstep,restore,full,itr,nem,&
                     gradt,eval,nref,dnrmg)
  !******************************************************************************!
  !                                                                              !
  !                           Update trust radius                                !
  !                                                                              !
  !******************************************************************************!
  implicit none

  ! External variables
  integer, intent(in   ) :: itr, nem
  integer, intent(  out) :: nref
  real*8,  intent(in   ) :: dnrmk, diffg, diffq, dnrmg
  real*8,  intent(inout) :: tradius
  logical, intent(in   ) :: full
  logical, intent(inout) :: rejected, trdef, nstep, restore

  real*8,  intent(in), pointer :: gradt(:), eval(:)

  ! Internal variables
  logical       :: flat
  logical, save :: was_flat = .false.
  character*6   :: cond
  integer       :: i, nred
  integer, save :: nrun  = 0, ncount(4) = 0

  real*8            :: R, tradold, dmaxa, dmaxn, pow
  real*8, save      :: diffg1, diffg2
  real*8, parameter :: d1q = 0.25d0, d3q = 0.75d0, d5q = 1.25d0

  !------------------------------------------------------------------------------!
  !                              Initialization                                  !
  !------------------------------------------------------------------------------!
  nred    = size(eval,1)
  tradold = tradius
  R       = diffg / diffq

  !------------------------------------------------------------------------------!
  !                            Trus-radius update                                !
  !------------------------------------------------------------------------------!
  if( dabs(diffg)<thresh .or. dabs(diffq)<thresh )then
    cond = '0     '
    if(rejected)then
      if(diffq > d0)then
        cond = '0r q>0'
        tradius  = tradius / 10.0d0
      else
        cond = '0r    '
        rejected = .false.
      endif
    endif
  else
    if(rejected)then
      if(diffq > d0)then
        cond = 'r  q>0'
        tradius  = tradius / 10.0d0
      else
        cond = 'r     '
        if(R < d0 ) cond = 'r  R<0'
        if(R < d1q) cond = 'r  R<1'
        if(R > d3q) cond = 'r  R>3'
        rejected = .false.
        if(R < d0 ) rejected = .true.
        if(R < d1q) tradius  = d1q* tradius
        if(R > d3q) tradius  = d2 * tradius
      endif
    else
      if(dnrmk > tradius)then
        if(full)then
          if(diffq < d0)then
            cond = 'Tq    '
            if(R < d0 ) cond = 'Tq R<0'
            if(R < d1q) cond = 'Tq R<1'
            if(R > d5q) cond = 'Tq R>5'
            tradius = dnrmk
            if(R < d0 ) rejected = .true.
            if(R < d1q) tradius  = d1q*tradius
            if(R > d5q) tradius  = d2*tradius
          else 
            cond = 'TQ    '
            if(R > d0 ) cond = 'TQ R>0'
            if((R < d1q).or.(R > d1+d3q)) cond = 'TQ R>7'
            if((R > d3q).and.(R < d5q))   cond = 'TQ R<5'
            tradius = dnrmk
            if(R > d0 ) rejected = .true.
            if((R < d1q).or.(R > d1+d3q)) tradius = d1q*tradius
            if((R > d3q).and.(R < d5q))   tradius = d2*tradius
          endif
        else
          cond = 'T  dia'
          rejected = .true.
        endif
      else
        if(diffq < d0)then
          cond = 'tq     '
          if(R < d0 ) cond = 'tq R<0'
          if(R < d1q) cond = 'tq R<1'
          if(R > d5q) cond = 'tq R>5'
          tradius = tradold
          if(trdef) tradius = dnrmk
          if(R < d0 ) rejected = .true.
          if(R < d1q) tradius  = d1q * dnrmk
          if(R > d5q) tradius  = max(d2*dnrmk,tradold)
        else 
          cond = 'tQ    '
          if(R > d0 ) cond = 'tQ R>0'
          if((R < d1q).or.(R > d1+d3q)) cond = 'tQ R>7'
          if((R > d3q).and.(R < d5q))   cond = 'tQ R<5'
          tradius = tradold
          if(trdef) tradius = dnrmk
          if(R > d0 ) rejected = .true.
          if((R < d1q).or.(R > d1+d3q)) tradius = d1q * dnrmk
          if((R > d3q).and.(R < d5q))   tradius = max(d2*dnrmk,tradold)
        endif
      endif
    endif
  endif

  !------------------------------------------------------------------------------!
  !                       Are we close to the minimum?                           !
  !            If both dmaxa, dmaxn are small we are near minimum.               !
  !------------------------------------------------------------------------------!
  dmaxn = d0
  do i=1,nred
    if(eval(i) > d0) exit
    dmaxn = dmaxn + dabs(d2*eval(i)*dnrmk)
  enddo
  if(nem == 0)then
    dmaxn = d1
  else
    dmaxn = dmaxn / dfloat(nem)
  endif
  dmaxa = dnrmg*dnrmk
  if(mnf%iprint.ge.2)then
    write(mnf%l,'(a,i4,2d14.6)') ' itr,dmaxa,dmaxn: ',itr,dmaxa,dmaxn
  endif

  !------------------------------------------------------------------------------!
  !         If all gradients in the direction of negative eigenvalues            !
  !              are below threshold, choose the best step from                  !
  !         (the purpose of this condition is to replace ineffective             !
  !                   rejection of Newton step - get_mu )                        !
  ! Newton step,                                        -->  rejected = false    !
  ! step on the Trust-region                            -->  rejected = true     !
  ! and step in the direction of positive eigenvalues.  -->  nstep    = false    !
  !------------------------------------------------------------------------------!
  nref = 0
  do i=1,nred
    if(dabs(gradt(i)) > thresh)then
      nref = i
      exit
    endif
  enddo
  if(nref == 0)then
    if(nem == 0)then
      write(mnf%l,'(a)')" WARNING: We are in minimum. Why I'm here??? "
      write(mnf%l,'(a)')" WARNING:    ==>  STOP "
      stop
    else
      nstep = .false.
    endif
  else
    if((eval(nref) > d0).and.(nem/=0)) nstep = .false.
    if( dmaxa<10.0d0*thresh .or. dmaxn<10.0d0*thresh ) nstep = .true.
  endif

  if( .not.nstep )then
    cond(3:3) = '!'
    if(restore)then          ! Newton step was restored
      rejected = .false.
      nstep    = .true.
      restore  = .false.
    elseif(nrun == 0)then    ! try side step
      diffg1   = diffg
      tradius  = tradold
      rejected = .true.
      nstep    = .false.
      nrun     = 1
    elseif( nrun == 1 .and. nref == 0 )then
      if( diffg1<diffg )then  ! Newton step                           
        tradius  = tradold                                                                    
        rejected = .false.                                                                    
        nstep    = .true.                                                                     
        restore  = .true.                                                                     
        nrun     = 0                                                                          
      else                    ! side step                             
        tradius  = tradold                                                                    
        rejected = .true.                                                                     
        nstep    = .false.                                                                    
        nrun     = 3                                                                          
      endif
    elseif(nrun == 1)then    ! try step on the Trust-region
      diffg2   = diffg
      tradius  = tradold
      rejected = .true.
      nstep    = .true.
      nrun     = 2
    elseif(nrun == 2)then    ! choose the best step
      if( diffg<d0 .or. diffg1<d0 .or. diffg2<d0 )then
        if( diffg<diffg1 .and. diffg<diffg2 )then       ! step on the Trust-region
          rejected = .false.
          nstep    = .true.
          nrun     = 0
        elseif( diffg1<diffg .and. diffg1<diffg2 )then  ! Newton step
          tradius  = tradold
          rejected = .false.
          nstep    = .true.
          restore  = .true.
          nrun     = 0
        elseif( diffg2<diffg .and. diffg2<diffg1 )then  ! side step
          tradius  = tradold
          rejected = .true.
          nstep    = .false.
          nrun     = 3
        endif
      else
        write(mnf%l,'(a)')" WARNING: Both regular Newton step and side step were rejected! "
        write(mnf%l,'(a)')" WARNING:    ==>  STOP "
        stop
      endif
    elseif(nrun == 3)then
      rejected = .false.
      nstep    = .true.
      nrun     = 0
    endif
  endif

  !------------------------------------------------------------------------------!
  !                       Check if we are in the flat area                       !
  !                                                                              !
  !          reject the step if diffg and diffq are both negative and            !
  !               lower than 1.0d-11 three and more times in row                 !
  !               lower than 1.0d-10 four  and more times in row                 !
  !               lower than 1.0d-09 five  and more times in row                 !
  !               lower than 1.0d-08 six   and more times in row                 !
  !------------------------------------------------------------------------------!
  flat = .false.
  if( diffg<d0 .and. diffq<d0 .and. .not.rejected .and. cond(3:3) /= '!' )then
    pow = 100.d0
    do i=1,size(ncount,1)
      if( dabs(diffg)<pow*thresh .and. dabs(diffq)<pow*thresh )then
        ncount(i) = ncount(i) + 1
      else
        ncount(i) = 0
      endif
      pow = pow * 10.d0
      if( ncount(i) > i+1 ) flat = .true.
    enddo
  else
    ncount = 0
  endif

  ! reject the step if flat and enlarge radius if we are still in the flat area
  if( flat )then
    if(was_flat) tradius = 10.0d0*tradius
    rejected  = .true.
    was_flat  = .true.
    cond(3:3) = '.'
  else
    was_flat  = .false.
  endif

  !------------------------------------------------------------------------------!
  !                Check if we are in the saddle point or maximum                !
  !------------------------------------------------------------------------------!
  if(dmaxa<thresh .and. nem/=0 .and. cond(3:3) /= '!')then
    rejected  = .true.
    if(cond(3:3) == ' ') cond(3:3) = 's'
    if(cond(3:3) == '.') cond(3:3) = 'd'
  endif

  if(tradius /= tradold) trdef = .false.

  write(mnf%l,'(/a,2x,a4, 4x,a, 6x,a, 5x,a, 4x,a, 3x,a, 4x,a, 5x,a, 6x,a)') &
    '  ','ITR','T_old','T_new ','|KAPPA|','  R   ','REJ','CON','DIFFG','DIFFQ'
  write(mnf%l,'(a,2x,i4,3d11.2,f11.5,3x,l1,3x,a,4d11.2)') &
    '*T',itr+1,tradold,tradius,dnrmk,R,rejected,cond,diffg,diffq


  end subroutine updateT

end module second_order_minimization
