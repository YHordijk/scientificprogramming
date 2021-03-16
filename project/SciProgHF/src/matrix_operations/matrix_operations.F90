
! radovan.bast@uit.no (2009-06-03):
! this module is the interface between matrix_genop_old and DIRAC
! please mind that the naming is common with Dalton
! so some things make only sense in the larger picture (Dalton and DIRAC)

! radovan: please do not use (also in the f90 terms) this module
!          directly, use this functionality only via module matrix_defop_old
!          otherwise you risk memory corruption!

module matrix_operations

  use matrix_stat
  use matrix_module

  implicit none

  public matrix
  public mat_init_old
  public mat_free
  public mat_mul
  public mat_zero
  public mat_scal
  public mat_trans
  public mat_dotproduct
  public mat_trab
  public mat_daxpy
  public mat_assign
  public mat_to_full
  public mat_set_from_full
  public mat_add
  public mat_copy
  public get_common_nz

  private

#include "dcbbas.h"
#include "dcborb.h"

contains

   function get_common_nz()
      integer :: get_common_nz
#include "dgroup.h"
      get_common_nz = nz
   end function

  subroutine mat_copy(alpha,a,b)
     implicit none
     REAL(8),  INTENT(IN)    :: alpha
     TYPE(Matrix), INTENT(IN)    :: a
     TYPE(Matrix), INTENT(INOUT) :: b
     integer                     :: i

!    do i = 1,a%nrow*a%ncol
!      b%elms(i) = alpha*a%elms(i)
!    enddo

     call mat_assign(b,a)
     if (alpha.ne.1.0d0) call mat_scal(alpha,b)

  end subroutine

  subroutine mat_add(alpha,a,beta,b,c)
     implicit none
     TYPE(Matrix), intent(IN) :: a, b
     REAL(8), INTENT(IN)  :: alpha, beta
     TYPE(Matrix)             :: c
     integer                  :: i

!    do i = 1,a%nrow*a%ncol
!      c%elms(i) = alpha*a%elms(i) + beta*b%elms(i)
!    enddo

     call mat_copy(alpha,a,c)
     call mat_daxpy(beta,b,c)

  end subroutine

  subroutine mat_set_from_full(afull,alpha,a, mat_label,unres3)
     implicit none
     real(8), INTENT(IN) :: afull(*)
     real(8), intent(in) :: alpha
     TYPE(Matrix)            :: a
     integer                 :: i
         character(*), INTENT(IN), OPTIONAL :: mat_label
         logical,intent(in) , OPTIONAL :: unres3

!    do i = 1,a%nrow*a%ncol
!      a%elms(i) = alpha*afull(i)
!    enddo

     i = a%nrow*a%ncol
     call dcopy (i,afull,1,a%elms,1)
     if (alpha.ne.1.0d0) call mat_scal(alpha,a)

  end subroutine

  subroutine mat_to_full(a, alpha, afull, mat_label)


     implicit none
     TYPE(Matrix), intent(in):: a
     real(8), intent(in) :: alpha
     real(8), intent(out):: afull(*)
     integer                 :: i
         character(*), INTENT(IN), OPTIONAL :: mat_label

!    do i = 1,a%nrow*a%ncol
!      afull(i) = alpha*a%elms(i)
!    enddo

     i = a%nrow*a%ncol*a%nz
     call dcopy (i,a%elms,1,afull,1)
     if (alpha.ne.1.0d0) call dscal(i,alpha,afull,1)

  end subroutine


  subroutine mat_assign(a,b)
     implicit none
     TYPE(Matrix), INTENT(INOUT) :: a
     TYPE(Matrix), INTENT(IN)    :: b
     integer                     :: i

!    do i = 1,a%nrow*a%ncol
!      a%elms(i) = b%elms(i)
!    enddo

     i = a%nrow*a%ncol*a%nz
     call dcopy (i,b%elms,1,a%elms,1)

  end subroutine

  subroutine mat_init_old(A, nrow, ncol, nz)

!   ----------------------------------------------------------------------------
    type(matrix)        :: A
    integer,           intent(in) :: nrow, ncol
    integer, optional, intent(in) :: nz
!   ----------------------------------------------------------------------------

    nullify(A%elms)
    nullify(A%iaux)
    nullify(A%raux)

    A%nrow = nrow
    A%ncol = ncol
    if(present(nz))then
      A%nz = nz
    else
      A%nz = get_common_nz()
    endif

    allocate(A%elms(A%nrow*A%ncol*A%nz*merge(2, 1, A%complex)))
    call add_allocation()

  end subroutine

  subroutine mat_free(A)

!   ----------------------------------------------------------------------------
    type(matrix) :: A
!   ----------------------------------------------------------------------------

    if (.not. associated(A%elms)) then
      call quit('error in mat_free: memory previously released')
    end if

    deallocate(A%elms)
    nullify(A%elms)
    call remove_allocation()

    if (associated(A%iaux)) deallocate(A%iaux)
    if (associated(A%raux)) deallocate(A%raux)

  end subroutine

  subroutine mat_zero(A)

!   ----------------------------------------------------------------------------
    type(matrix), intent(inout) :: A
!   ----------------------------------------------------------------------------

    A%elms = 0.0d0

  end subroutine

  subroutine mat_mul(A, B, transa, transb, alpha, beta, C)

!   C      = alpha*A*B + beta*C
!   C_{mn} = ... A_{mk} B_{kn} ...
!   transa/transb: 'T' or 't' (transposed), 'N' or 'n' (normal)
!   ----------------------------------------------------------------------------
    type(matrix), intent(in)    :: A, B
    character,    intent(in)    :: transa, transb
    real(8),      intent(in)    :: alpha, beta
    type(matrix), intent(inout) :: C
!   ----------------------------------------------------------------------------
    integer                     :: m, n, k, l
    character                   :: ta, tb
    integer                     :: iq(4) = (/1, 2, 3, 4/)
!   ----------------------------------------------------------------------------
#include "dgroup.h"

    if (nbsym > 1) then
      call iqpack(ipqtoq(1, A%irep), A%nz, &
                  ipqtoq(1, B%irep), B%nz, &
                  iq               , C%nz)
    end if

    select case (transa)
      case ('n', 'N')
        m  = A%nrow
        k  = A%ncol
        ta = 'N'
      case ('t', 'T')
        m  = A%ncol
        k  = A%nrow
        ta = 'H'
      case default
        call quit('unknown trans possibility in mat_mul')
    end select

    select case (transb)
      case ('n', 'N')
        l  = B%nrow
        n  = B%ncol
        tb = 'N'
      case ('t', 'T')
        l  = B%ncol
        n  = B%nrow
        tb = 'H'
      case default
        call quit('unknown trans possibility in mat_mul')
    end select

    if (k /= l) then
      call quit('dimensions do not match in mat_mul')
    end if

    call qgemm(m, n, k,                                                  &
               alpha,                                                    &
               ta, 'N', ipqtoq(1, A%irep), A%elms, A%nrow, A%ncol, A%nz, &
               tb, 'N', ipqtoq(1, B%irep), B%elms, B%nrow, B%ncol, B%nz, &
               beta,    iq               , C%elms, C%nrow, C%ncol, C%nz)

  end subroutine

  subroutine mat_scal(alpha, A)

!   ----------------------------------------------------------------------------
    real(8),      intent(in)    :: alpha
    type(matrix), intent(inout) :: A
!   ----------------------------------------------------------------------------

    call dscal(A%nrow*A%ncol*A%nz, alpha, A%elms, 1)

  end subroutine

  subroutine mat_trans(A, B)

!   ----------------------------------------------------------------------------
    type(matrix), intent(in) :: A
    type(matrix)             :: B
!   ----------------------------------------------------------------------------
    integer                  :: i, j, i_dim, j_dim, iz
    real(8)                  :: q
!   ----------------------------------------------------------------------------
#include "dgroup.h"

    if (B%nrow /= A%ncol .or. B%ncol /= A%nrow) then
      call quit('wrong dimensions in mat_trans')
    end if

    i_dim = A%nrow
    j_dim = A%ncol

    iz = 1
    do iz = 1, A%nz

      q = 1.0d0
      if (ipqtoq(iz, A%irep) > 1) q = -1.0d0

      do j = 1, j_dim
        do i = 1, i_dim
              B%elms((iz - 1)*i_dim*j_dim + i_dim*(i - 1) + j) &
          = q*A%elms((iz - 1)*i_dim*j_dim + i_dim*(j - 1) + i)
        end do
      end do
    end do

  end subroutine

  function mat_dotproduct(A, B) result(r)

!   ----------------------------------------------------------------------------
    type(matrix), intent(in) :: A, B
    real(8)                  :: r
!   ----------------------------------------------------------------------------
    real(8),      external   :: ddot
!   ----------------------------------------------------------------------------

    if (A%nrow*A%ncol /= B%nrow*B%ncol) then
      call quit('wrong dimensions in mat_dotproduct')
    end if

    r = 0.0d0

    if (ieor(A%irep, B%irep) == 0) then
      r = ddot(A%nrow*A%ncol*A%nz, A%elms, 1, B%elms, 1)
    end if

  end function

  function mat_trab(A, B) result(r)

!   ----------------------------------------------------------------------------
    type(matrix), intent(in) :: A, B
    real(8)                  :: r
!   ----------------------------------------------------------------------------
    real(8),      external   :: ddot
    real(8)                  :: q
    integer                  :: i_dim, j_dim, j, iz
!   ----------------------------------------------------------------------------
#include "dgroup.h"

    if (A%ncol /= B%nrow .or. A%nrow /= B%ncol) then
      call quit('wrong dimensions in mat_trab')
    end if

    r = 0.0d0

    i_dim = A%nrow
    j_dim = A%ncol

    if (ieor(A%irep, B%irep) == 0) then
       do iz = 1, B%nz

         q = 1.0d0
         if (ipqtoq(iz, B%irep) > 1) q = -1.0d0

         do j = 1, j_dim
           r = r + q*ddot(i_dim, A%elms((iz - 1)*i_dim*j_dim + i_dim*(j - 1) + 1), 1, &
                                 B%elms((iz - 1)*i_dim*j_dim + j), j_dim)
         enddo
       end do
    end if

  end function

  subroutine mat_daxpy(alpha, X, Y)

!   Y = alpha*X + Y
!   ----------------------------------------------------------------------------
    real(8),      intent(in)    :: alpha
    type(matrix), intent(in)    :: X
    type(matrix), intent(inout) :: Y
!   ----------------------------------------------------------------------------

    if (X%nrow /= Y%nrow .or. X%ncol /= Y%ncol) then
      call quit('wrong dimensions in mat_daxpy')
    endif

    call daxpy(X%nrow*X%ncol*X%nz, alpha, X%elms, 1, Y%elms,1)

  end subroutine

end module
