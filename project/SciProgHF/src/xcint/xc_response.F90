module xc_response

#if defined MOD_OPENRSP
  use xc_derv
  use density_eval
  use openrsp_cfg
  use xcfun_1_0
  use xc_max_block_length

  implicit none

  public xc_response_contribution

  private
    real(8), parameter     :: tiny = 1.0d-30
    integer, parameter     :: max_nr_dmat = 20

!    "the universe works on a math equation
!     that never even ever really ends in the end"
!                      modest mouse - never ending math equation - building nothing out of something


!     we have five density variables:
!       n
!       s
!       z = \nabla n \cdot nabla n
!       y = \nabla n \cdot nabla s
!       x = \nabla s \cdot nabla s
!
!     response contribution blocks:
!         n    s    z    y    x
!        _________________________
!     n | n    ns   nz   rest rest
!     s | ns   ns   sz   rest rest
!     z | nz   sz   nz   rest rest
!     y | rest rest rest rest rest
!     x | rest rest rest rest rest
!
!     \nabla n contribution | s contribution | \nabla s contribution | blocks
!     -----------------------------------------------------------------------
!     no                    | no             | no                    | n
!     no                    | yes            | no                    | n + ns
!     yes                   | no             | no                    | n + nz
!     yes                   | yes            | no                    | n + nz + ns + sz
!     yes                   | yes            | yes                   | n + nz + ns + sz + rest

  logical :: calculate_blocks_n    = .true.
  logical :: calculate_blocks_ns   = .false.
  logical :: calculate_blocks_nz   = .true.
  logical :: calculate_blocks_sz   = .false.
  logical :: calculate_blocks_rest = .false.

contains

  subroutine xc_response_contribution(order,     &
                                      mat_dim,   &
                                      nz,        &
                                      fmat,      &
                                      nr_dmat_p, &
                                      dmat_p,    &
                                      isym_d,    &
                                      ih_d,    &
                                      isym_f,    &
                                      use_gga,   &
                                      ao, &
                                      n_0,       &
                                      gn_0,      &
                                      px,        &
                                      py,        &
                                      pz,        &
                                      w,         &
                                      buffer,      &
                                      derv_length, &
                                      derv)

!   ----------------------------------------------------------------------------
    integer, intent(in)    :: order
    integer, intent(in)    :: mat_dim
    integer, intent(in)    :: nz
    real(8), intent(inout) :: fmat(mat_dim, mat_dim, nz, *)
    integer, intent(in)    :: nr_dmat_p
    real(8), intent(in)    :: dmat_p(mat_dim, mat_dim, nz, nr_dmat_p)
    integer, intent(in)    :: isym_d(nr_dmat_p)
    integer, intent(in)    :: ih_d(nr_dmat_p)
    integer, intent(in)    :: isym_f(*)
    logical, intent(in)    :: use_gga
    real(8), intent(in)    :: ao(*)
    real(8), intent(in)    :: n_0
    real(8), intent(in)    :: gn_0(3)
    real(8), intent(in)    :: px
    real(8), intent(in)    :: py
    real(8), intent(in)    :: pz
    real(8), intent(in)    :: w
    real(8), intent(in)    :: buffer(*)
    integer, intent(in)    :: derv_length
    real(8)                :: derv(derv_length, max_block_length)
!   ----------------------------------------------------------------------------
    real(8)                :: u
    real(8)                :: v(3)

    real(8)                :: n(0:max_nr_dmat)
    real(8)                :: gn(3, 0:max_nr_dmat)
    real(8)                :: z(0:max_nr_dmat, 0:max_nr_dmat)

    integer                :: i, j
    logical                :: u_above_threshold
    logical                :: v_above_threshold

    real(8), external      :: ddot
    real(8)                :: xcfun_in(5, max_block_length)
!   ----------------------------------------------------------------------------

#ifdef INT_STAR8
    if (use_gga) then
       print *, 'controlled stop: only int32'
       print *, 'index integer type conversion fails with int64 for GGA'
       stop 1
    end if
#endif

!   get perturbed densities
!   get perturbed densitiy gradients and gradient norms

    n(0) = n_0
    gn(:, 0) = gn_0
    do i = 1, nr_dmat_p
       if (use_gga) then
          if (ih_d(i) == 0) then
             call get_gn_nonhermitian(n(i),  &
                         gn(1, i),           &
                         isym_d(i) - 1,      &
                         mat_dim,            &
                         dmat_p(1, 1, 1, i), &
                         buffer,             &
                         ao)
          else
             call get_gn(n(i),               &
                         gn(1, i),           &
                         isym_d(i) - 1,      &
                         mat_dim,            &
                         dmat_p(1, 1, 1, i), &
                         buffer,             &
                         ao)
          end if
       else
          call get_n(n(i),               &
                     isym_d(i) - 1,      &
                     mat_dim,            &
                     dmat_p(1, 1, 1, i), &
                     buffer,             &
                     ao)
       end if
    end do

!   otherwise undefined with (a)lda
    z(0, 0) = 0.0d0

    if (use_gga) then
      do i = 0, nr_dmat_p
        do j = i, nr_dmat_p
!         the perturbed z carries factor 2
!         because derivative of z gives 2 identical terms
!         (one "left", one "right")
          z(i, j) = 2.0d0*ddot(3, gn(1, i), 1, gn(1, j), 1)
        end do
      end do
!     correct wrong prefactor 2 of the unperturbed z
      z(0, 0) = 0.5d0*z(0, 0)

    end if

    xcfun_in = 0.0d0
    xcfun_in(1, 1) = n(0)
    xcfun_in(3, 1) = z(0, 0)

    derv = 0.0d0
    call xc_eval(xc_fun%id, order, 1, xcfun_in, derv)

    u = 0.0d0
    v = 0.0d0

    select case (order)

      case (2)
        call xc_response_2(u, v, derv, n, gn, z, use_gga)

      case (3)
        call xc_response_3(u, v, derv, n, gn, z, use_gga)

      case (4)
        call xc_response_4(u, v, derv, n, gn, z, use_gga)

      case (5)
        call xc_response_5(u, v, derv, n, gn, z, use_gga)

      case default
        call quit('xc_response_contribution: order not implemented')

    end select

!   multiply by the integration weight
    u = u*w
    v = v*w

    u_above_threshold = .false.
    if (dabs(u) > tiny) then
      u_above_threshold = .true.
    end if

    v_above_threshold = .false.
    if (maxval((/dabs(v(1)), dabs(v(2)), dabs(v(3))/)) > tiny) then
      v_above_threshold = .true.
    end if

    if (v_above_threshold) then
       call nabla_omega_real(fmat,          &
                             1,             &
                             isym_f(1) - 1, &
                             mat_dim,       &
                             nz,            &
                             u,             &
                             2.0d0*v,       &
                             ao)
    else
       if (u_above_threshold) then
          call omega_real(fmat,          &
                          1,             &
                          isym_f(1) - 1, &
                          mat_dim,       &
                          nz,            &
                          u,             &
                          ao)
       end if
    end if

  end subroutine

  subroutine xc_response_2(u, v, derv, n, gn, z, use_gga)

!   ----------------------------------------------------------------------------
    real(8), intent(inout) :: u
    real(8), intent(inout) :: v(3)
    real(8), intent(in)    :: derv(*)
    real(8), intent(in)    :: n(0:max_nr_dmat)
    real(8), intent(in)    :: gn(3, 0:max_nr_dmat)
    real(8), intent(in)    :: z(0:max_nr_dmat, 0:max_nr_dmat)
    logical, intent(in)    :: use_gga
!   ----------------------------------------------------------------------------
    real(8)                :: tu
    real(8)                :: tv(3)
!   ----------------------------------------------------------------------------

    tu = n(1)
    u  = u + derv(xc_index(xc_fun%id, (/2, 0, 0, 0, 0/)))*tu
    if (use_gga) then
    v  = v + derv(xc_index(xc_fun%id, (/1, 0, 1, 0, 0/)))*tu*gn(:, 0)
    tu = z(0, 1)
    v  = v + derv(xc_index(xc_fun%id, (/0, 0, 2, 0, 0/)))*tu*gn(:, 0)
    u  = u + derv(xc_index(xc_fun%id, (/1, 0, 1, 0, 0/)))*tu
    tv = gn(:, 1)
    v  = v + derv(xc_index(xc_fun%id, (/0, 0, 1, 0, 0/)))*tv(:)
    end if

  end subroutine

  subroutine xc_response_3(u, v, derv, n, gn, z, use_gga)

!   ----------------------------------------------------------------------------
    real(8), intent(inout) :: u
    real(8), intent(inout) :: v(3)
    real(8), intent(in)    :: derv(*)
    real(8), intent(in)    :: n(0:max_nr_dmat)
    real(8), intent(in)    :: gn(3, 0:max_nr_dmat)
    real(8), intent(in)    :: z(0:max_nr_dmat, 0:max_nr_dmat)
    logical, intent(in)    :: use_gga
!   ----------------------------------------------------------------------------
    real(8)                :: tu
    real(8)                :: tv(3)
!   ----------------------------------------------------------------------------

    tu = n(1)*n(2)
    u  = u + derv(xc_index(xc_fun%id, (/3, 0, 0, 0, 0/)))*tu
    if (use_gga) then
    v  = v + derv(xc_index(xc_fun%id, (/2, 0, 1, 0, 0/)))*tu*gn(:, 0)
    tv = z(0, 1)*gn(:, 2) &
       + gn(:, 1)*z(0, 2)
    v  = v + derv(xc_index(xc_fun%id, (/0, 0, 2, 0, 0/)))*tv(:)
    tu = n(1)*z(0, 2) &
       + z(0, 1)*n(2)
    v  = v + derv(xc_index(xc_fun%id, (/1, 0, 2, 0, 0/)))*tu*gn(:, 0)
    u  = u + derv(xc_index(xc_fun%id, (/2, 0, 1, 0, 0/)))*tu
    tu = z(0, 1)*z(0, 2)
    v  = v + derv(xc_index(xc_fun%id, (/0, 0, 3, 0, 0/)))*tu*gn(:, 0)
    u  = u + derv(xc_index(xc_fun%id, (/1, 0, 2, 0, 0/)))*tu
    tv = n(1)*gn(:, 2) &
       + gn(:, 1)*n(2)
    v  = v + derv(xc_index(xc_fun%id, (/1, 0, 1, 0, 0/)))*tv(:)
    tu = z(1, 2)
    v  = v + derv(xc_index(xc_fun%id, (/0, 0, 2, 0, 0/)))*tu*gn(:, 0)
    u  = u + derv(xc_index(xc_fun%id, (/1, 0, 1, 0, 0/)))*tu
    end if

  end subroutine

  subroutine xc_response_4(u, v, derv, n, gn, z, use_gga)

!   ----------------------------------------------------------------------------
    real(8), intent(inout) :: u
    real(8), intent(inout) :: v(3)
    real(8), intent(in)    :: derv(*)
    real(8), intent(in)    :: n(0:max_nr_dmat)
    real(8), intent(in)    :: gn(3, 0:max_nr_dmat)
    real(8), intent(in)    :: z(0:max_nr_dmat, 0:max_nr_dmat)
    logical, intent(in)    :: use_gga
!   ----------------------------------------------------------------------------
    real(8)                :: tu
    real(8)                :: tv(3)
!   ----------------------------------------------------------------------------

    tu = n(1)*n(2)*n(3)
    u  = u + derv(xc_index(xc_fun%id, (/4, 0, 0, 0, 0/)))*tu
    if (use_gga) then
    v  = v + derv(xc_index(xc_fun%id, (/3, 0, 1, 0, 0/)))*tu*gn(:, 0)
    end if
    tu = n(4)*n(3) &
       + n(5)*n(2) &
       + n(6)*n(1)
    u  = u + derv(xc_index(xc_fun%id, (/3, 0, 0, 0, 0/)))*tu
    if (use_gga) then
    v  = v + derv(xc_index(xc_fun%id, (/2, 0, 1, 0, 0/)))*tu*gn(:, 0)
    tv = z(0, 4)*gn(:, 3) &
       + z(1, 2)*gn(:, 3) &
       + z(0, 5)*gn(:, 2) &
       + z(1, 3)*gn(:, 2) &
       + z(0, 6)*gn(:, 1) &
       + z(2, 3)*gn(:, 1) &
       + gn(:, 4)*z(0, 3) &
       + gn(:, 5)*z(0, 2) &
       + gn(:, 6)*z(0, 1)
    v  = v + derv(xc_index(xc_fun%id, (/0, 0, 2, 0, 0/)))*tv(:)
    tu = n(1)*n(2)*z(0, 3) &
       + n(1)*z(0, 2)*n(3) &
       + z(0, 1)*n(2)*n(3)
    v  = v + derv(xc_index(xc_fun%id, (/2, 0, 2, 0, 0/)))*tu*gn(:, 0)
    u  = u + derv(xc_index(xc_fun%id, (/3, 0, 1, 0, 0/)))*tu
    tu = n(4)*z(0, 3) &
       + n(5)*z(0, 2) &
       + n(6)*z(0, 1) &
       + z(0, 4)*n(3) &
       + z(1, 2)*n(3) &
       + z(0, 5)*n(2) &
       + z(1, 3)*n(2) &
       + z(0, 6)*n(1) &
       + z(2, 3)*n(1)
    v  = v + derv(xc_index(xc_fun%id, (/1, 0, 2, 0, 0/)))*tu*gn(:, 0)
    u  = u + derv(xc_index(xc_fun%id, (/2, 0, 1, 0, 0/)))*tu
    tu = z(0, 1)*z(0, 2)*z(0, 3)
    v  = v + derv(xc_index(xc_fun%id, (/0, 0, 4, 0, 0/)))*tu*gn(:, 0)
    u  = u + derv(xc_index(xc_fun%id, (/1, 0, 3, 0, 0/)))*tu
    tu = z(0, 4)*z(0, 3) &
       + z(1, 2)*z(0, 3) &
       + z(0, 5)*z(0, 2) &
       + z(1, 3)*z(0, 2) &
       + z(0, 6)*z(0, 1) &
       + z(2, 3)*z(0, 1)
    v  = v + derv(xc_index(xc_fun%id, (/0, 0, 3, 0, 0/)))*tu*gn(:, 0)
    u  = u + derv(xc_index(xc_fun%id, (/1, 0, 2, 0, 0/)))*tu
    tv = n(4)*gn(:, 3) &
       + n(5)*gn(:, 2) &
       + n(6)*gn(:, 1) &
       + gn(:, 4)*n(3) &
       + gn(:, 5)*n(2) &
       + gn(:, 6)*n(1)
    v  = v + derv(xc_index(xc_fun%id, (/1, 0, 1, 0, 0/)))*tv(:)
    tv = n(1)*z(0, 2)*gn(:, 3) &
       + n(1)*gn(:, 2)*z(0, 3) &
       + z(0, 1)*n(2)*gn(:, 3) &
       + z(0, 1)*gn(:, 2)*n(3) &
       + gn(:, 1)*n(2)*z(0, 3) &
       + gn(:, 1)*z(0, 2)*n(3)
    v  = v + derv(xc_index(xc_fun%id, (/1, 0, 2, 0, 0/)))*tv(:)
    tu = n(1)*z(0, 2)*z(0, 3) &
       + z(0, 1)*n(2)*z(0, 3) &
       + z(0, 1)*z(0, 2)*n(3)
    v  = v + derv(xc_index(xc_fun%id, (/1, 0, 3, 0, 0/)))*tu*gn(:, 0)
    u  = u + derv(xc_index(xc_fun%id, (/2, 0, 2, 0, 0/)))*tu
    tv = n(1)*n(2)*gn(:, 3) &
       + n(1)*gn(:, 2)*n(3) &
       + gn(:, 1)*n(2)*n(3)
    v  = v + derv(xc_index(xc_fun%id, (/2, 0, 1, 0, 0/)))*tv(:)
    tu = z(1, 6) &
       + z(2, 5) &
       + z(3, 4)
    v  = v + derv(xc_index(xc_fun%id, (/0, 0, 2, 0, 0/)))*tu*gn(:, 0)
    u  = u + derv(xc_index(xc_fun%id, (/1, 0, 1, 0, 0/)))*tu
    tv = z(0, 1)*z(0, 2)*gn(:, 3) &
       + z(0, 1)*gn(:, 2)*z(0, 3) &
       + gn(:, 1)*z(0, 2)*z(0, 3)
    v  = v + derv(xc_index(xc_fun%id, (/0, 0, 3, 0, 0/)))*tv(:)
    end if

  end subroutine

  subroutine xc_response_5(u, v, derv, n, gn, z, use_gga)

!   ----------------------------------------------------------------------------
    real(8), intent(inout) :: u
    real(8), intent(inout) :: v(3)
    real(8), intent(in)    :: derv(*)
    real(8), intent(in)    :: n(0:max_nr_dmat)
    real(8), intent(in)    :: gn(3, 0:max_nr_dmat)
    real(8), intent(in)    :: z(0:max_nr_dmat, 0:max_nr_dmat)
    logical, intent(in)    :: use_gga
!   ----------------------------------------------------------------------------
    real(8)                :: tu
    real(8)                :: tv(3)
!   ----------------------------------------------------------------------------

    tu = n(1)*n(2)*n(3)*n(4)
    u  = u + derv(xc_index(xc_fun%id, (/5, 0, 0, 0, 0/)))*tu
    if (use_gga) then
    v  = v + derv(xc_index(xc_fun%id, (/4, 0, 1, 0, 0/)))*tu*gn(:, 0)
    end if
    tu = n(11)*n(4) &
       + n(12)*n(3) &
       + n(13)*n(2) &
       + n(14)*n(1) &
       + n(5)*n(10) &
       + n(6)*n(9) &
       + n(7)*n(8)
    u  = u + derv(xc_index(xc_fun%id, (/3, 0, 0, 0, 0/)))*tu
    if (use_gga) then
    v  = v + derv(xc_index(xc_fun%id, (/2, 0, 1, 0, 0/)))*tu*gn(:, 0)
    end if
    tu = n(5)*n(3)*n(4) &
       + n(6)*n(2)*n(4) &
       + n(7)*n(2)*n(3) &
       + n(8)*n(1)*n(4) &
       + n(9)*n(1)*n(3) &
       + n(10)*n(1)*n(2)
    u  = u + derv(xc_index(xc_fun%id, (/4, 0, 0, 0, 0/)))*tu
    if (use_gga) then
    v  = v + derv(xc_index(xc_fun%id, (/3, 0, 1, 0, 0/)))*tu*gn(:, 0)
    tv = z(0, 11)*gn(:, 4) &
       + z(1, 8)*gn(:, 4) &
       + z(2, 6)*gn(:, 4) &
       + z(3, 5)*gn(:, 4) &
       + z(0, 12)*gn(:, 3) &
       + z(1, 9)*gn(:, 3) &
       + z(2, 7)*gn(:, 3) &
       + z(4, 5)*gn(:, 3) &
       + z(0, 13)*gn(:, 2) &
       + z(1, 10)*gn(:, 2) &
       + z(3, 7)*gn(:, 2) &
       + z(4, 6)*gn(:, 2) &
       + z(0, 14)*gn(:, 1) &
       + z(2, 10)*gn(:, 1) &
       + z(3, 9)*gn(:, 1) &
       + z(4, 8)*gn(:, 1) &
       + gn(:, 11)*z(0, 4) &
       + gn(:, 12)*z(0, 3) &
       + gn(:, 13)*z(0, 2) &
       + gn(:, 14)*z(0, 1) &
       + z(0, 5)*gn(:, 10) &
       + z(1, 2)*gn(:, 10) &
       + z(0, 6)*gn(:, 9) &
       + z(1, 3)*gn(:, 9) &
       + z(0, 7)*gn(:, 8) &
       + z(1, 4)*gn(:, 8) &
       + z(0, 8)*gn(:, 7) &
       + z(2, 3)*gn(:, 7) &
       + z(0, 9)*gn(:, 6) &
       + z(2, 4)*gn(:, 6) &
       + z(0, 10)*gn(:, 5) &
       + z(3, 4)*gn(:, 5)
    v  = v + derv(xc_index(xc_fun%id, (/0, 0, 2, 0, 0/)))*tv(:)
    tv = n(11)*gn(:, 4) &
       + n(12)*gn(:, 3) &
       + n(13)*gn(:, 2) &
       + n(14)*gn(:, 1) &
       + gn(:, 11)*n(4) &
       + gn(:, 12)*n(3) &
       + gn(:, 13)*n(2) &
       + gn(:, 14)*n(1) &
       + n(5)*gn(:, 10) &
       + n(6)*gn(:, 9) &
       + n(7)*gn(:, 8) &
       + n(8)*gn(:, 7) &
       + n(9)*gn(:, 6) &
       + n(10)*gn(:, 5)
    v  = v + derv(xc_index(xc_fun%id, (/1, 0, 1, 0, 0/)))*tv(:)
    tu = n(1)*z(0, 2)*z(0, 3)*z(0, 4) &
       + z(0, 1)*n(2)*z(0, 3)*z(0, 4) &
       + z(0, 1)*z(0, 2)*n(3)*z(0, 4) &
       + z(0, 1)*z(0, 2)*z(0, 3)*n(4)
    v  = v + derv(xc_index(xc_fun%id, (/1, 0, 4, 0, 0/)))*tu*gn(:, 0)
    u  = u + derv(xc_index(xc_fun%id, (/2, 0, 3, 0, 0/)))*tu
    tv = z(0, 5)*z(0, 3)*gn(:, 4) &
       + z(1, 2)*z(0, 3)*gn(:, 4) &
       + z(0, 6)*z(0, 2)*gn(:, 4) &
       + z(1, 3)*z(0, 2)*gn(:, 4) &
       + z(0, 7)*z(0, 2)*gn(:, 3) &
       + z(1, 4)*z(0, 2)*gn(:, 3) &
       + z(0, 8)*z(0, 1)*gn(:, 4) &
       + z(2, 3)*z(0, 1)*gn(:, 4) &
       + z(0, 9)*z(0, 1)*gn(:, 3) &
       + z(2, 4)*z(0, 1)*gn(:, 3) &
       + z(0, 10)*z(0, 1)*gn(:, 2) &
       + z(3, 4)*z(0, 1)*gn(:, 2) &
       + z(0, 5)*gn(:, 3)*z(0, 4) &
       + z(1, 2)*gn(:, 3)*z(0, 4) &
       + z(0, 6)*gn(:, 2)*z(0, 4) &
       + z(1, 3)*gn(:, 2)*z(0, 4) &
       + z(0, 7)*gn(:, 2)*z(0, 3) &
       + z(1, 4)*gn(:, 2)*z(0, 3) &
       + z(0, 8)*gn(:, 1)*z(0, 4) &
       + z(2, 3)*gn(:, 1)*z(0, 4) &
       + z(0, 9)*gn(:, 1)*z(0, 3) &
       + z(2, 4)*gn(:, 1)*z(0, 3) &
       + z(0, 10)*gn(:, 1)*z(0, 2) &
       + z(3, 4)*gn(:, 1)*z(0, 2) &
       + gn(:, 5)*z(0, 3)*z(0, 4) &
       + gn(:, 6)*z(0, 2)*z(0, 4) &
       + gn(:, 7)*z(0, 2)*z(0, 3) &
       + gn(:, 8)*z(0, 1)*z(0, 4) &
       + gn(:, 9)*z(0, 1)*z(0, 3) &
       + gn(:, 10)*z(0, 1)*z(0, 2)
    v  = v + derv(xc_index(xc_fun%id, (/0, 0, 3, 0, 0/)))*tv(:)
    tu = z(0, 1)*z(0, 2)*z(0, 3)*z(0, 4)
    v  = v + derv(xc_index(xc_fun%id, (/0, 0, 5, 0, 0/)))*tu*gn(:, 0)
    u  = u + derv(xc_index(xc_fun%id, (/1, 0, 4, 0, 0/)))*tu
    tu = z(0, 11)*z(0, 4) &
       + z(1, 8)*z(0, 4) &
       + z(2, 6)*z(0, 4) &
       + z(3, 5)*z(0, 4) &
       + z(0, 12)*z(0, 3) &
       + z(1, 9)*z(0, 3) &
       + z(2, 7)*z(0, 3) &
       + z(4, 5)*z(0, 3) &
       + z(0, 13)*z(0, 2) &
       + z(1, 10)*z(0, 2) &
       + z(3, 7)*z(0, 2) &
       + z(4, 6)*z(0, 2) &
       + z(0, 14)*z(0, 1) &
       + z(2, 10)*z(0, 1) &
       + z(3, 9)*z(0, 1) &
       + z(4, 8)*z(0, 1) &
       + z(0, 5)*z(0, 10) &
       + z(0, 5)*z(3, 4) &
       + z(1, 2)*z(0, 10) &
       + z(1, 2)*z(3, 4) &
       + z(0, 6)*z(0, 9) &
       + z(0, 6)*z(2, 4) &
       + z(1, 3)*z(0, 9) &
       + z(1, 3)*z(2, 4) &
       + z(0, 7)*z(0, 8) &
       + z(0, 7)*z(2, 3) &
       + z(1, 4)*z(0, 8) &
       + z(1, 4)*z(2, 3)
    v  = v + derv(xc_index(xc_fun%id, (/0, 0, 3, 0, 0/)))*tu*gn(:, 0)
    u  = u + derv(xc_index(xc_fun%id, (/1, 0, 2, 0, 0/)))*tu
    tv = n(5)*n(3)*gn(:, 4) &
       + n(6)*n(2)*gn(:, 4) &
       + n(7)*n(2)*gn(:, 3) &
       + n(8)*n(1)*gn(:, 4) &
       + n(9)*n(1)*gn(:, 3) &
       + n(10)*n(1)*gn(:, 2) &
       + n(5)*gn(:, 3)*n(4) &
       + n(6)*gn(:, 2)*n(4) &
       + n(7)*gn(:, 2)*n(3) &
       + n(8)*gn(:, 1)*n(4) &
       + n(9)*gn(:, 1)*n(3) &
       + n(10)*gn(:, 1)*n(2) &
       + gn(:, 5)*n(3)*n(4) &
       + gn(:, 6)*n(2)*n(4) &
       + gn(:, 7)*n(2)*n(3) &
       + gn(:, 8)*n(1)*n(4) &
       + gn(:, 9)*n(1)*n(3) &
       + gn(:, 10)*n(1)*n(2)
    v  = v + derv(xc_index(xc_fun%id, (/2, 0, 1, 0, 0/)))*tv(:)
    tu = n(5)*z(0, 3)*z(0, 4) &
       + n(6)*z(0, 2)*z(0, 4) &
       + n(7)*z(0, 2)*z(0, 3) &
       + n(8)*z(0, 1)*z(0, 4) &
       + n(9)*z(0, 1)*z(0, 3) &
       + n(10)*z(0, 1)*z(0, 2) &
       + z(0, 5)*n(3)*z(0, 4) &
       + z(1, 2)*n(3)*z(0, 4) &
       + z(0, 6)*n(2)*z(0, 4) &
       + z(1, 3)*n(2)*z(0, 4) &
       + z(0, 7)*n(2)*z(0, 3) &
       + z(1, 4)*n(2)*z(0, 3) &
       + z(0, 8)*n(1)*z(0, 4) &
       + z(2, 3)*n(1)*z(0, 4) &
       + z(0, 9)*n(1)*z(0, 3) &
       + z(2, 4)*n(1)*z(0, 3) &
       + z(0, 10)*n(1)*z(0, 2) &
       + z(3, 4)*n(1)*z(0, 2) &
       + z(0, 5)*z(0, 3)*n(4) &
       + z(1, 2)*z(0, 3)*n(4) &
       + z(0, 6)*z(0, 2)*n(4) &
       + z(1, 3)*z(0, 2)*n(4) &
       + z(0, 7)*z(0, 2)*n(3) &
       + z(1, 4)*z(0, 2)*n(3) &
       + z(0, 8)*z(0, 1)*n(4) &
       + z(2, 3)*z(0, 1)*n(4) &
       + z(0, 9)*z(0, 1)*n(3) &
       + z(2, 4)*z(0, 1)*n(3) &
       + z(0, 10)*z(0, 1)*n(2) &
       + z(3, 4)*z(0, 1)*n(2)
    v  = v + derv(xc_index(xc_fun%id, (/1, 0, 3, 0, 0/)))*tu*gn(:, 0)
    u  = u + derv(xc_index(xc_fun%id, (/2, 0, 2, 0, 0/)))*tu
    tv = n(1)*z(0, 2)*z(0, 3)*gn(:, 4) &
       + n(1)*z(0, 2)*gn(:, 3)*z(0, 4) &
       + n(1)*gn(:, 2)*z(0, 3)*z(0, 4) &
       + z(0, 1)*n(2)*z(0, 3)*gn(:, 4) &
       + z(0, 1)*n(2)*gn(:, 3)*z(0, 4) &
       + z(0, 1)*z(0, 2)*n(3)*gn(:, 4) &
       + z(0, 1)*z(0, 2)*gn(:, 3)*n(4) &
       + z(0, 1)*gn(:, 2)*n(3)*z(0, 4) &
       + z(0, 1)*gn(:, 2)*z(0, 3)*n(4) &
       + gn(:, 1)*n(2)*z(0, 3)*z(0, 4) &
       + gn(:, 1)*z(0, 2)*n(3)*z(0, 4) &
       + gn(:, 1)*z(0, 2)*z(0, 3)*n(4)
    v  = v + derv(xc_index(xc_fun%id, (/1, 0, 3, 0, 0/)))*tv(:)
    tu = z(0, 5)*z(0, 3)*z(0, 4) &
       + z(1, 2)*z(0, 3)*z(0, 4) &
       + z(0, 6)*z(0, 2)*z(0, 4) &
       + z(1, 3)*z(0, 2)*z(0, 4) &
       + z(0, 7)*z(0, 2)*z(0, 3) &
       + z(1, 4)*z(0, 2)*z(0, 3) &
       + z(0, 8)*z(0, 1)*z(0, 4) &
       + z(2, 3)*z(0, 1)*z(0, 4) &
       + z(0, 9)*z(0, 1)*z(0, 3) &
       + z(2, 4)*z(0, 1)*z(0, 3) &
       + z(0, 10)*z(0, 1)*z(0, 2) &
       + z(3, 4)*z(0, 1)*z(0, 2)
    v  = v + derv(xc_index(xc_fun%id, (/0, 0, 4, 0, 0/)))*tu*gn(:, 0)
    u  = u + derv(xc_index(xc_fun%id, (/1, 0, 3, 0, 0/)))*tu
    tu = n(1)*n(2)*n(3)*z(0, 4) &
       + n(1)*n(2)*z(0, 3)*n(4) &
       + n(1)*z(0, 2)*n(3)*n(4) &
       + z(0, 1)*n(2)*n(3)*n(4)
    v  = v + derv(xc_index(xc_fun%id, (/3, 0, 2, 0, 0/)))*tu*gn(:, 0)
    u  = u + derv(xc_index(xc_fun%id, (/4, 0, 1, 0, 0/)))*tu
    tv = z(0, 1)*z(0, 2)*z(0, 3)*gn(:, 4) &
       + z(0, 1)*z(0, 2)*gn(:, 3)*z(0, 4) &
       + z(0, 1)*gn(:, 2)*z(0, 3)*z(0, 4) &
       + gn(:, 1)*z(0, 2)*z(0, 3)*z(0, 4)
    v  = v + derv(xc_index(xc_fun%id, (/0, 0, 4, 0, 0/)))*tv(:)
    tv = n(5)*z(0, 3)*gn(:, 4) &
       + n(6)*z(0, 2)*gn(:, 4) &
       + n(7)*z(0, 2)*gn(:, 3) &
       + n(8)*z(0, 1)*gn(:, 4) &
       + n(9)*z(0, 1)*gn(:, 3) &
       + n(10)*z(0, 1)*gn(:, 2) &
       + n(5)*gn(:, 3)*z(0, 4) &
       + n(6)*gn(:, 2)*z(0, 4) &
       + n(7)*gn(:, 2)*z(0, 3) &
       + n(8)*gn(:, 1)*z(0, 4) &
       + n(9)*gn(:, 1)*z(0, 3) &
       + n(10)*gn(:, 1)*z(0, 2) &
       + z(0, 5)*n(3)*gn(:, 4) &
       + z(1, 2)*n(3)*gn(:, 4) &
       + z(0, 6)*n(2)*gn(:, 4) &
       + z(1, 3)*n(2)*gn(:, 4) &
       + z(0, 7)*n(2)*gn(:, 3) &
       + z(1, 4)*n(2)*gn(:, 3) &
       + z(0, 8)*n(1)*gn(:, 4) &
       + z(2, 3)*n(1)*gn(:, 4) &
       + z(0, 9)*n(1)*gn(:, 3) &
       + z(2, 4)*n(1)*gn(:, 3) &
       + z(0, 10)*n(1)*gn(:, 2) &
       + z(3, 4)*n(1)*gn(:, 2) &
       + z(0, 5)*gn(:, 3)*n(4) &
       + z(1, 2)*gn(:, 3)*n(4) &
       + z(0, 6)*gn(:, 2)*n(4) &
       + z(1, 3)*gn(:, 2)*n(4) &
       + z(0, 7)*gn(:, 2)*n(3) &
       + z(1, 4)*gn(:, 2)*n(3) &
       + z(0, 8)*gn(:, 1)*n(4) &
       + z(2, 3)*gn(:, 1)*n(4) &
       + z(0, 9)*gn(:, 1)*n(3) &
       + z(2, 4)*gn(:, 1)*n(3) &
       + z(0, 10)*gn(:, 1)*n(2) &
       + z(3, 4)*gn(:, 1)*n(2) &
       + gn(:, 5)*n(3)*z(0, 4) &
       + gn(:, 6)*n(2)*z(0, 4) &
       + gn(:, 7)*n(2)*z(0, 3) &
       + gn(:, 8)*n(1)*z(0, 4) &
       + gn(:, 9)*n(1)*z(0, 3) &
       + gn(:, 10)*n(1)*z(0, 2) &
       + gn(:, 5)*z(0, 3)*n(4) &
       + gn(:, 6)*z(0, 2)*n(4) &
       + gn(:, 7)*z(0, 2)*n(3) &
       + gn(:, 8)*z(0, 1)*n(4) &
       + gn(:, 9)*z(0, 1)*n(3) &
       + gn(:, 10)*z(0, 1)*n(2)
    v  = v + derv(xc_index(xc_fun%id, (/1, 0, 2, 0, 0/)))*tv(:)
    tu = n(11)*z(0, 4) &
       + n(12)*z(0, 3) &
       + n(13)*z(0, 2) &
       + n(14)*z(0, 1) &
       + z(0, 11)*n(4) &
       + z(1, 8)*n(4) &
       + z(2, 6)*n(4) &
       + z(3, 5)*n(4) &
       + z(0, 12)*n(3) &
       + z(1, 9)*n(3) &
       + z(2, 7)*n(3) &
       + z(4, 5)*n(3) &
       + z(0, 13)*n(2) &
       + z(1, 10)*n(2) &
       + z(3, 7)*n(2) &
       + z(4, 6)*n(2) &
       + z(0, 14)*n(1) &
       + z(2, 10)*n(1) &
       + z(3, 9)*n(1) &
       + z(4, 8)*n(1) &
       + n(5)*z(0, 10) &
       + n(5)*z(3, 4) &
       + n(6)*z(0, 9) &
       + n(6)*z(2, 4) &
       + n(7)*z(0, 8) &
       + n(7)*z(2, 3) &
       + n(8)*z(0, 7) &
       + n(8)*z(1, 4) &
       + n(9)*z(0, 6) &
       + n(9)*z(1, 3) &
       + n(10)*z(0, 5) &
       + n(10)*z(1, 2)
    v  = v + derv(xc_index(xc_fun%id, (/1, 0, 2, 0, 0/)))*tu*gn(:, 0)
    u  = u + derv(xc_index(xc_fun%id, (/2, 0, 1, 0, 0/)))*tu
    tu = n(5)*n(3)*z(0, 4) &
       + n(6)*n(2)*z(0, 4) &
       + n(7)*n(2)*z(0, 3) &
       + n(8)*n(1)*z(0, 4) &
       + n(9)*n(1)*z(0, 3) &
       + n(10)*n(1)*z(0, 2) &
       + n(5)*z(0, 3)*n(4) &
       + n(6)*z(0, 2)*n(4) &
       + n(7)*z(0, 2)*n(3) &
       + n(8)*z(0, 1)*n(4) &
       + n(9)*z(0, 1)*n(3) &
       + n(10)*z(0, 1)*n(2) &
       + z(0, 5)*n(3)*n(4) &
       + z(1, 2)*n(3)*n(4) &
       + z(0, 6)*n(2)*n(4) &
       + z(1, 3)*n(2)*n(4) &
       + z(0, 7)*n(2)*n(3) &
       + z(1, 4)*n(2)*n(3) &
       + z(0, 8)*n(1)*n(4) &
       + z(2, 3)*n(1)*n(4) &
       + z(0, 9)*n(1)*n(3) &
       + z(2, 4)*n(1)*n(3) &
       + z(0, 10)*n(1)*n(2) &
       + z(3, 4)*n(1)*n(2)
    v  = v + derv(xc_index(xc_fun%id, (/2, 0, 2, 0, 0/)))*tu*gn(:, 0)
    u  = u + derv(xc_index(xc_fun%id, (/3, 0, 1, 0, 0/)))*tu
    tv = n(1)*n(2)*z(0, 3)*gn(:, 4) &
       + n(1)*n(2)*gn(:, 3)*z(0, 4) &
       + n(1)*z(0, 2)*n(3)*gn(:, 4) &
       + n(1)*z(0, 2)*gn(:, 3)*n(4) &
       + n(1)*gn(:, 2)*n(3)*z(0, 4) &
       + n(1)*gn(:, 2)*z(0, 3)*n(4) &
       + z(0, 1)*n(2)*n(3)*gn(:, 4) &
       + z(0, 1)*n(2)*gn(:, 3)*n(4) &
       + z(0, 1)*gn(:, 2)*n(3)*n(4) &
       + gn(:, 1)*n(2)*n(3)*z(0, 4) &
       + gn(:, 1)*n(2)*z(0, 3)*n(4) &
       + gn(:, 1)*z(0, 2)*n(3)*n(4)
    v  = v + derv(xc_index(xc_fun%id, (/2, 0, 2, 0, 0/)))*tv(:)
    tu = n(1)*n(2)*z(0, 3)*z(0, 4) &
       + n(1)*z(0, 2)*n(3)*z(0, 4) &
       + n(1)*z(0, 2)*z(0, 3)*n(4) &
       + z(0, 1)*n(2)*n(3)*z(0, 4) &
       + z(0, 1)*n(2)*z(0, 3)*n(4) &
       + z(0, 1)*z(0, 2)*n(3)*n(4)
    v  = v + derv(xc_index(xc_fun%id, (/2, 0, 3, 0, 0/)))*tu*gn(:, 0)
    u  = u + derv(xc_index(xc_fun%id, (/3, 0, 2, 0, 0/)))*tu
    tv = n(1)*n(2)*n(3)*gn(:, 4) &
       + n(1)*n(2)*gn(:, 3)*n(4) &
       + n(1)*gn(:, 2)*n(3)*n(4) &
       + gn(:, 1)*n(2)*n(3)*n(4)
    v  = v + derv(xc_index(xc_fun%id, (/3, 0, 1, 0, 0/)))*tv(:)
    tu = z(1, 14) &
       + z(2, 13) &
       + z(3, 12) &
       + z(4, 11) &
       + z(5, 10) &
       + z(6, 9) &
       + z(7, 8)
    v  = v + derv(xc_index(xc_fun%id, (/0, 0, 2, 0, 0/)))*tu*gn(:, 0)
    u  = u + derv(xc_index(xc_fun%id, (/1, 0, 1, 0, 0/)))*tu
    end if

  end subroutine
#endif

end module
