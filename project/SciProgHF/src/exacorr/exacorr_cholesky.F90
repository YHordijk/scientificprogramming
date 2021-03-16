module exacorr_cholesky

!routines to perform cholesky decomposition

!       written by Johann Pototschnig, Summer 2019
!       some refactoring, Lucas Visscher, april 2020

use, intrinsic:: ISO_C_BINDING

implicit none

complex(8), parameter :: ZERO=(0.D0,0.D0),ONE=(1.D0,0.D0),MINUS_ONE=(-1.D0,0.D0)

private

!talsh routines:
public talsh_get_diag
public talsh_find_shells
public talsh_get_shell
public talsh_find_el
public talsh_square_el
public talsh_set_mask_done

contains

subroutine talsh_get_diag(Diag)

! routine to compute diagonal 2D tensor

  use tensor_algebra
  use talsh
  use exacorr_datatypes
  use exacorr_global
  use module_interest_eri

  implicit none

  type(talsh_tens_t), intent(inout) :: Diag

  real(8), pointer     :: DmR(:,:)
  complex(8),pointer   :: DmC(:,:)
  type(C_PTR)          :: body_p
  integer              :: rank
  integer              :: dims2(1:2)
  integer              :: ierr
  integer              :: DataKind(1),nData

  type(basis_func_info_t), allocatable :: gto(:)  ! arrays with basis function information
  integer                              :: nao, nshells
  integer                              :: basis_angular

  integer :: lq, lp, nq, np, q, p
  integer :: i, j, ijij
  integer :: ioff, joff
  real(8) :: eq,ep,cq,cp
  real(8) :: xq,yq,zq,xp,yp,zp
  real(8),save :: gout(21*21*21*21) ! hardwired for maximum  l-value of 5 ? (s=1,p=3,d=6,f=10,g=15,h=21)
  real(8), parameter :: FIJKL = 1.0
  
  ! get acces to the tensor after determining its type and dimensions
  rank=talsh_tensor_rank(Diag)
  if (rank.ne.2) stop 'error: wrong rank in find_shells_talsh'
  ierr = talsh_tensor_dimensions(Diag,rank,dims2)
  if (ierr.ne.0) stop 'error: wrong dimension in find_shells_talsh'
  ierr=talsh_tensor_data_kind(Diag,nData,DataKind)
  if (ierr.ne.0) stop 'error in getting DataKind'
  ierr=talsh_tensor_get_body_access(Diag,body_p,DataKind(1),0,DEV_HOST)
  select case (DataKind(1))
  case (R8)
     call c_f_pointer(body_p,DmR,dims2)
  case (C8)
     call c_f_pointer(body_p,DmC,dims2)
  case default
     print*, "wrong datakind :",DataKind(1)
  end select

  ! get basis information
  nao = get_nao()
  call get_gtos(1,nao,gto,nshells)
  basis_angular = get_basis_angular()

  ! compute diagonal elements
  joff = 0
  do q = 1, nshells
     lq   =  gto(q)%orb_momentum
     eq   =  gto(q)%exponent(1)
     xq   =  gto(q)%coord(1)
     yq   =  gto(q)%coord(2)
     zq   =  gto(q)%coord(3)
     cq   =  gto(q)%coefficient(1)
     nq   =  nfunctions(lq,basis_angular)

     ioff = 0
     do p = 1, nshells
        lp   =  gto(p)%orb_momentum
        ep   =  gto(p)%exponent(1)
        xp   =  gto(p)%coord(1)
        yp   =  gto(p)%coord(2)
        zp   =  gto(p)%coord(3)
        cp   =  gto(p)%coefficient(1)
        np   =  nfunctions(lp,basis_angular)

        !get eri values
        call interest_eri('llll',FIJKL,gout,&
                          lp,ep,xp,yp,zp,cp,&
                          lq,eq,xq,yq,zq,cq,&
                          lp,ep,xp,yp,zp,cp,&
                          lq,eq,xq,yq,zq,cq )

        do j = 1, nq
          do i = 1, np
            ijij = (j-1)*np*nq*np+(i-1)*nq*np+(j-1)*np+i
            if (DataKind(1)==R8) then
               DmR(i+ioff,j+joff) = gout(ijij)
            else
               DmC(i+ioff,j+joff) = dcmplx(gout(ijij),0.D0)
            end if
          end do
        end do

        ioff = ioff + np
     end do
     joff = joff + nq
  end do

end subroutine talsh_get_diag

subroutine talsh_find_shells(Diag, Filter, lshell, kshell, D_max)

! routine to find shell with largest element

  use tensor_algebra
  use talsh
  use exacorr_datatypes
  use exacorr_global

  implicit none

  type(talsh_tens_t), intent(inout) :: Diag
  type(talsh_tens_t), intent(inout) :: Filter
  integer, intent(out)              :: lshell, kshell
  complex(8), intent(out)           :: D_max


  integer              :: nData,DataKind(1)
  real(8), pointer     :: DmR(:,:)
  complex(8), pointer  :: DmC(:,:)
  real(4), pointer     :: rfilt(:,:)
  logical, allocatable :: lfilt(:,:)
  type(C_PTR)          :: body_p
  integer              :: rank
  integer              :: dims2(1:2), i_max(1:2)
  complex(8)           :: D_tot
  integer              :: ierr

  rank=talsh_tensor_rank(Diag)
  if (rank.ne.2) stop 'error: wrong rank in find_shells_talsh'
  ierr = talsh_tensor_dimensions(Diag,rank,dims2)
  if (ierr.ne.0) stop 'error: wrong dimension in find_shells_talsh'
  ierr=talsh_tensor_data_kind(Diag,nData,DataKind)
  if (ierr.ne.0) stop 'error in getting DataKind'
  ierr=talsh_tensor_get_body_access(Diag,body_p,DataKind(1),0,DEV_HOST)
  select case (DataKind(1))
  case (R8)
     call c_f_pointer(body_p,DmR,dims2)
  case (C8)
     call c_f_pointer(body_p,DmC,dims2)
  case default
     print*, "wrong datakind :",DataKind(1)
  end select

  rank=talsh_tensor_rank(Filter)
  if (rank.ne.2) stop 'error: wrong rank in find_shells_talsh'
  ierr = talsh_tensor_dimensions(Filter,rank,dims2)
  if (ierr.ne.0) stop 'error: wrong dimension in find_shells_talsh'
  ierr=talsh_tensor_get_body_access(Filter,body_p,R4,0,DEV_HOST)
  call c_f_pointer(body_p,rfilt,dims2)
  allocate(lfilt(dims2(1),dims2(2)))
  lfilt=rfilt.gt.0

  if (DataKind(1)==R8) then
     !find largest element
     i_max=maxloc(abs(DmR))
     D_tot=DmR(i_max(1), i_max(2))

     !find largest non treated element
     i_max=maxloc(abs(DmR), mask=lfilt)
     D_max=dcmplx(DmR(i_max(1), i_max(2)),0.D0)
  else
     !find largest element
     i_max=maxloc(abs(DmC))
     D_tot=DmC(i_max(1), i_max(2))

     !find largest non treated element
     i_max=maxloc(abs(DmC), mask=lfilt)
     D_max=DmC(i_max(1), i_max(2))
  end if

  if (abs(D_tot).gt.abs(D_max)) then
    print *, 'Warning: error in cholesky will not be smaller than ', D_tot
  end if

  lshell=get_shell_index(i_max(1))
  kshell=get_shell_index(i_max(2))

  deallocate(lfilt)

end subroutine talsh_find_shells

subroutine talsh_get_shell(Dshell, lsh, ksh)

! routine to compute all ao integrals for a shell pair

  use tensor_algebra
  use talsh
  use exacorr_datatypes
  use exacorr_global
  use module_interest_eri

  implicit none 

  integer, intent(in)            :: lsh, ksh
  type(talsh_tens_t), intent(inout) :: Dshell

  real(8), pointer     :: DlkR(:,:,:,:)
  complex(8), pointer  :: DlkC(:,:,:,:)
  type(C_PTR)          :: body_p
  integer              :: rank
  integer              :: dims4(1:4)
  integer              :: ierr


  type(basis_func_info_t), allocatable :: gto(:)  ! arrays with basis function information
  integer                              :: nao, nshells
  integer                              :: basis_angular
  integer                              :: nData,DataKind(1)

  !variables required to compute 2el integrals
  integer(8)         :: ijkl
  integer(8)         :: ioff,joff,koff,loff
  integer(8)         :: ish,jsh
  integer(8)         :: kshell,lshell
  integer            :: i, ni, li
  real(8)            :: ei, ci, xi, yi, zi
  integer            :: j, nj, lj
  real(8)            :: ej, cj, xj, yj, zj
  integer            :: k, nk, lk
  real(8)            :: ek, ck, xk, yk, zk
  integer            :: l, nl, ll
  real(8)            :: el, cl, xl, yl, zl
  real(8), parameter :: FIJKL = 1.0
  integer, parameter :: LMAX = 21 ! hardwired for maximum  l-value of 5 ? (s=1,p=3,d=6,f=10,g=15,h=21)
  real(8), save      :: gout(LMAX*LMAX*LMAX*LMAX) 
  real(8)            :: limit

  ! get basis information
  nao = get_nao()
  call get_gtos(1,nao,gto,nshells)
  basis_angular = get_basis_angular()

  !set parameters for the fixed shells
  loff=0
  do lshell = 1, lsh-1
    ll   =  gto(lshell)%orb_momentum
    nl   =  nfunctions(ll,basis_angular)
    loff = loff + nl
  end do
  ll   =  gto(lsh)%orb_momentum
  el   =  gto(lsh)%exponent(1)
  xl   =  gto(lsh)%coord(1)
  yl   =  gto(lsh)%coord(2)
  zl   =  gto(lsh)%coord(3)
  cl   =  gto(lsh)%coefficient(1)
  nl   =  nfunctions(ll,basis_angular)

  koff=0
  do kshell = 1, ksh-1
    lk   =  gto(kshell)%orb_momentum
    nk   =  nfunctions(lk,basis_angular)
    koff = koff + nk
  end do
  lk   =  gto(ksh)%orb_momentum
  ek   =  gto(ksh)%exponent(1)
  xk   =  gto(ksh)%coord(1)
  yk   =  gto(ksh)%coord(2)
  zk   =  gto(ksh)%coord(3)
  ck   =  gto(ksh)%coefficient(1)
  nk   =  nfunctions(lk,basis_angular)

  ! get acces to the tensor
  rank=talsh_tensor_rank(Dshell)
  if (rank.ne.4) stop 'error: wrong rank in get_shell_talsh'
  ierr = talsh_tensor_dimensions(Dshell,rank,dims4)
  if (ierr.ne.0) stop 'error: wrong dimension in get_shell_talsh'
  ierr=talsh_tensor_data_kind(Dshell,nData,DataKind)
  if (ierr.ne.0) stop 'error in getting DataKind'
  ierr=talsh_tensor_get_body_access(Dshell,body_p,DataKind(1),0,DEV_HOST)
  select case (DataKind(1))
  case (R8)
     call c_f_pointer(body_p,DlkR,dims4)
  case (C8)
     call c_f_pointer(body_p,DlkC,dims4)
  case default
     print*, "wrong datakind :",DataKind(1)
  end select

  !loop over all the shells for the other indices
  joff = 0
  do jsh = 1, nshells
    lj   =  gto(jsh)%orb_momentum
    ej   =  gto(jsh)%exponent(1)
    xj   =  gto(jsh)%coord(1)
    yj   =  gto(jsh)%coord(2)
    zj   =  gto(jsh)%coord(3)
    cj   =  gto(jsh)%coefficient(1)
    nj   =  nfunctions(lj,basis_angular)

    ioff = 0
    do ish = 1, nshells
      li   =  gto(ish)%orb_momentum
      ei   =  gto(ish)%exponent(1)
      xi   =  gto(ish)%coord(1)
      yi   =  gto(ish)%coord(2)
      zi   =  gto(ish)%coord(3)
      ci   =  gto(ish)%coefficient(1)
      ni   =  nfunctions(li,basis_angular)

!     output order of eri is (5,c,d,5,a,b), so input k,l first and then i,j to get the order that we want
      call interest_eri('llll',FIJKL,gout,&
                        lk,ek,xk,yk,zk,ck,&
                        ll,el,xl,yl,zl,cl,&
                        li,ei,xi,yi,zi,ci,&
                        lj,ej,xj,yj,zj,cj )

      !save all elements for this shell pair
      do k = 1, nk
        do l = 1, nl
          do j = 1, nj
            do i = 1, ni
              ijkl = (l-1)*nk*nj*ni+(k-1)*nj*ni+(j-1)*ni+i
              if (DataKind(1)==R8) then
                 DlkR(i+ioff,j+joff,l,k) = gout(ijkl)
              else
                 DlkC(i+ioff,j+joff,l,k) = dcmplx(gout(ijkl),0.D0)
              end if
            end do
          end do
        end do
      end do

      ioff = ioff + ni
    end do

    joff = joff + nj
  end do

end subroutine talsh_get_shell

subroutine talsh_find_el(Diag, Filter, lsh, ksh, l, k, D_el)

! routine to locate largest element in a shell

  use tensor_algebra
  use talsh
  use exacorr_datatypes
  use exacorr_global

  implicit none
  
  integer, intent(in)               :: lsh, ksh
  type(talsh_tens_t), intent(inout) :: Filter
  type(talsh_tens_t), intent(inout) :: Diag
  integer, intent(out)              :: l, k
  complex(8), intent(out)           :: D_el

  real(8), pointer     :: DmR(:,:)
  complex(8), pointer  :: DmC(:,:)
  type(C_PTR)          :: body_p
  integer              :: rank
  integer              :: dims2(1:2), mind(1:2)
  integer              :: ierr
  real(4), pointer     :: rfilt(:,:)
  logical, allocatable :: lfilt(:,:)
  integer              :: nData,DataKind(1)

  type(basis_func_info_t), allocatable :: gto(:)  ! arrays with basis function information
  integer                              :: nao, nshells
  integer                              :: basis_angular

  integer(8)         :: koff, loff
  integer            :: kshell, lshell
  integer            :: nk, nl
  integer            :: ll, lk

  logical, allocatable :: Dmask(:,:)

  ! get basis information
  nao = get_nao()
  call get_gtos(1,nao,gto,nshells)
  basis_angular = get_basis_angular()

  !get shell information
  loff=0
  do lshell = 1, lsh-1
    ll   =  gto(lshell)%orb_momentum
    nl   =  nfunctions(ll,basis_angular)
    loff = loff + nl
  end do
  ll   =  gto(lsh)%orb_momentum
  nl   =  nfunctions(ll,basis_angular)

  koff=0
  do kshell = 1, ksh-1
    lk   =  gto(kshell)%orb_momentum
    nk   =  nfunctions(lk,basis_angular)
    koff = koff + nk
  end do
  lk   =  gto(ksh)%orb_momentum
  nk   =  nfunctions(lk,basis_angular)

  ! get acces to the tensor
  rank=talsh_tensor_rank(Diag)
  if (rank.ne.2) stop 'error: wrong rank in get_shell_talsh'
  ierr = talsh_tensor_dimensions(Diag,rank,dims2)
  if (ierr.ne.0) stop 'error: wrong dimension in get_shell_talsh'
  ierr=talsh_tensor_data_kind(Diag,nData,DataKind)
  if (ierr.ne.0) stop 'error in getting DataKind'
  ierr=talsh_tensor_get_body_access(Diag,body_p,DataKind(1),0,DEV_HOST)
  select case (DataKind(1))
  case (R8)
     call c_f_pointer(body_p,DmR,dims2)
  case (C8)
     call c_f_pointer(body_p,DmC,dims2)
  case default
     print*, "wrong datakind :",DataKind(1)
  end select

  ! set up Mask
  rank=talsh_tensor_rank(Filter)
  if (rank.ne.2) stop 'error: wrong rank in get_shell_talsh'
  ierr = talsh_tensor_dimensions(Filter,rank,dims2)
  if (ierr.ne.0) stop 'error: wrong dimension in get_shell_talsh'
  ierr=talsh_tensor_get_body_access(Filter,body_p,R4,0,DEV_HOST)
  call c_f_pointer(body_p,rfilt,dims2)
  allocate(lfilt(dims2(1),dims2(2)))
  lfilt=.false.

  allocate(Dmask(dims2(1),dims2(2)))
  Dmask=.false.
  Dmask(koff+1:koff+nk,loff+1:loff+nl)=.true.
  where (Dmask) lfilt = rfilt.gt.0

  if (any(lfilt)) then
    if (DataKind(1)==R8) then
       mind=maxloc(abs(DmR), mask=lfilt)
       k=mind(1)
       l=mind(2)
       D_el=dcmplx(DmR(l,k),0.D0)
     else
       mind=maxloc(abs(DmC), mask=lfilt)
       k=mind(1)
       l=mind(2)
       D_el=DmC(l,k)
     end if
  else
    k=0
    l=0
    D_el=ZERO
  end if
  
  deallocate(lfilt)
  deallocate(Dmask)

end subroutine talsh_find_el

subroutine talsh_square_el(h_tensor)

! squaring talsh tensor

  use tensor_algebra
  use talsh

  implicit none
  
  type(talsh_tens_t), intent(inout)    :: h_tensor

  integer  ::  dims2(1:2)
  integer  ::  dims3(1:3)
  integer  ::  dims4(1:4)
  integer  ::  ierr
  integer  ::  rank, rrank
  integer  ::  i, j, k, l
  type(C_PTR)         :: body_p
  complex(8), pointer :: tens2(:, :)
  complex(8), pointer :: tens3(:, :, :)
  complex(8), pointer :: tens4(:, :, :, :)

  rank=talsh_tensor_rank(h_tensor)
        
  if (rank.eq.2) then
    ierr = talsh_tensor_dimensions(h_tensor,rrank,dims2)
    if (ierr.ne.0 .or. rrank.ne.rank) stop 'error in square_el: wrong rank'

    ierr=talsh_tensor_get_body_access(h_tensor,body_p,C8,0,DEV_HOST)
    call c_f_pointer(body_p,tens2,dims2)
    
    do i = 1, dims2(1)
      do j = 1, dims2(2)
        tens2(i,j)=tens2(i,j)*tens2(i,j)
      end do
    end do

  else if (rank.eq.3) then
    ierr = talsh_tensor_dimensions(h_tensor,rrank,dims3)
    if (ierr.ne.0 .or. rrank.ne.rank) stop 'error in square_el: wrong rank'

    ierr=talsh_tensor_get_body_access(h_tensor,body_p,C8,0,DEV_HOST)
    call c_f_pointer(body_p,tens3,dims3)

    do i = 1, dims3(1)
      do j = 1, dims3(2)
        do k = 1, dims3(3)
          tens3(i, j, k)=tens3(i, j, k)*tens3(i, j, k)
        end do
      end do
    end do

  else if (rank.eq.4) then
    ierr = talsh_tensor_dimensions(h_tensor,rrank,dims4)
    if (ierr.ne.0 .or. rrank.ne.rank) stop 'error in square_el: wrong rank'

    ierr=talsh_tensor_get_body_access(h_tensor,body_p,C8,0,DEV_HOST)
    call c_f_pointer(body_p,tens4,dims4)

    do i = 1, dims4(1)
      do j = 1, dims4(2)
        do k = 1, dims4(3)
          do l = 1, dims4(4)
            tens4(i, j, k, l)=tens4(i, j, k, l)*tens4(i, j, k, l)
          end do
        end do
      end do
    end do

  else
    stop 'error in print_tensor: only ranks 2,3 and 4 implemented'
  end if

end subroutine talsh_square_el

subroutine talsh_set_mask_done(Filter,l,k)

! updating mask

  use tensor_algebra
  use talsh
  
  implicit none
  
  integer, intent(in)               :: l, k
  type(talsh_tens_t), intent(inout) :: Filter

  type(C_PTR)          :: body_p
  integer              :: dims2(1:2)
  integer              :: ierr
  real(4), pointer     :: rfilt(:,:)
  integer              :: rank

  ! set up Mask
  rank=talsh_tensor_rank(Filter)
  if (rank.ne.2) stop 'error: wrong rank in get_shell_talsh'
  ierr = talsh_tensor_dimensions(Filter,rank,dims2)
  if (ierr.ne.0) stop 'error: wrong dimension in get_shell_talsh'
  ierr=talsh_tensor_get_body_access(Filter,body_p,R4,0,DEV_HOST)
  call c_f_pointer(body_p,rfilt,dims2)
  
  rfilt(l,k)=-1.0

end subroutine talsh_set_mask_done

end module exacorr_cholesky
