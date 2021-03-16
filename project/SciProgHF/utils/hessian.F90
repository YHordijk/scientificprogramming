  program hessian
! This program reads information extracted from a GOSCI run and
! sets up the electronic Hessian
  integer :: i,j,ii,ia,ji,ja,stat
  integer :: indi,indj,indk,indl,nxop
  integer :: norb,nocc,nvrt,nvar,ntkr,nokr,nvkr
  integer, allocatable :: jxop(:,:),ind(:)
  real(8), allocatable :: orbeig(:)
  real(8), allocatable :: twoint(:,:,:,:)
  real(8), allocatable :: ablock(:,:),bblock(:,:)

! We start by getting the information extracted from GOSCI
! ========================================================
  open (5, Form='FORMATTED')
!
! First line contains total number of orbitals and occupied ones
! --------------------------------------------------------------
  read(5,*) norb,nocc
  nvrt=norb-nocc
! 
! Calculate the number of Kramers pair for each orbital class
! and make a mapping of indices
! Positive index is unbarred; negative index is barred
! -----------------------------------------------------------
  ntkr=norb/2
  nokr=nocc/2
  nvkr=nvrt/2
  allocate(ind(-ntkr:ntkr))
  do i=1,ntkr
    ind( i)= 2*i-1
    ind(-i)= 2*i
  enddo
!
! Next section is a list of orbital energies
! ------------------------------------------
  read(5,*)
  allocate (orbeig(ntkr))
  do i = 1,ntkr
    read(5,*) orbeig(i)
    read(5,*) 
  enddo
!  write(6,*) (orbeig(i),i=1,ntkr)
!
! Final section is the list of non-zero two electron integrals
! ------------------------------------------------------------
  allocate (twoint(norb,norb,norb,norb))
  twoint = 0
  read(5,*)
  stat = 0
  do while (stat.eq.0)
    read(5,*,iostat=stat) indi,indj,indk,indl,twoint(indi,indj,indk,indl)
  enddo
  close(5,status='KEEP')
!
! We set up orbital rotation indices 
! ========================================================================
  nvar=nocc*nvrt
  allocate (jxop(2,nvar))
  i = 0
!
! unbarred-unbarred
! -----------------------
  do ia = nokr+1,ntkr
    do ii = 1,nokr
      i = i+1
      jxop(1,i)=ind( ia)
      jxop(2,i)=ind( ii)
    enddo
  enddo
!
! unbarred-barred
! ---------------
  do ia = nokr+1,ntkr
    do ii = 1,nokr
      i = i+1
      jxop(1,i)=ind( ia)
      jxop(2,i)=ind(-ii)
    enddo
  enddo
!
! barred-barred
! ---------------
  do ia = nokr+1,ntkr
    do ii = 1,nokr
      i = i+1
      jxop(1,i)=ind(-ia)
      jxop(2,i)=ind(-ii)
    enddo
  enddo
!
! barred-unbarred
! ---------------
  do ia = nokr+1,ntkr
    do ii = 1,nokr
      i = i+1
      jxop(1,i)=ind(-ia)
      jxop(2,i)=ind( ii)
    enddo
  enddo
!
! Set up A-block of Hessian (coupling terms only)
! ===============================================
  allocate (ablock(nvar,nvar))
  do j=1,nvar
    ja=jxop(1,j)
    ji=jxop(2,j)
    do i=1,nvar
      ia=jxop(1,i)
      ii=jxop(2,i)
      ablock(i,j)=twoint(ia,ii,ji,ja)-twoint(ia,ja,ji,ii)
    enddo
  enddo
!  write(6,'(/A/)') '*** A block ***'
  do i = 1,nvar
    write(6,*) ablock(i,:)
  enddo
!
! Set up B-block of Hessian (coupling terms only)
! ===============================================
  allocate (bblock(nvar,nvar))
  do j=1,nvar
    ja=jxop(1,j)
    ji=jxop(2,j)
    do i=1,nvar
      ia=jxop(1,i)
      ii=jxop(2,i)
      bblock(i,j)=twoint(ia,ii,ja,ji)-twoint(ia,ji,ja,ii)
    enddo
  enddo
!  write(6,'(/A/)') '*** B block ***'
  do i = 1,nvar
    write(6,*) bblock(i,:)
  enddo

  end
