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

module xc_blas3

   use interface_ao
   use xc_max_block_length
   use xc_derv
!fixme dft_cfg should not be used here, move xc derv business to separate module
   use dft_cfg
   use xc_ac
   use extra

   implicit none

   public xc_contribution_blas3
   public get_block_threshold

   private

   integer, allocatable :: ao_compress_index(:)

contains

   subroutine xc_contribution_blas3(                         &
                                    is_gga,                  &
                                    is_gga_qr,               &
                                    is_gga_qi,               &
                                    do_lr,                   &
                                    block_length,            &
                                    block_threshold,         &
                                    mat_dim,                 &
                                    nz,                      &
                                    buffer,                  &
                                    ao_block,                &
                                    nr_fmat,                 &
                                    fmat,                    &
                                    dmat_0,                  &
                                    dmat_0_saop_outer,       &
                                    dmat_0_saop_resp,        &
                                    nr_dmat,                 &
                                    dmat_p,                  &
                                    dmat_pg_sym,             &
                                    dmat_ih_sym,             &
                                    derv_length,             &
                                    derv,                    &
                                    max_fun_derv,            &
                                    w,                       &
                                    nr_electrons_integrated, &
                                    xc_energy                &
                                   )

!     --------------------------------------------------------------------------
      logical, intent(in)    :: is_gga
      logical, intent(in)    :: is_gga_qr
      logical, intent(in)    :: is_gga_qi
      logical, intent(in)    :: do_lr
      integer, intent(in)    :: block_length
      real(8), intent(in)    :: block_threshold
      integer, intent(in)    :: mat_dim
      integer, intent(in)    :: nz
      real(8)                :: buffer(*)
      real(8), intent(in)    :: ao_block(block_length, *)
      integer, intent(in)    :: nr_fmat
      real(8), intent(inout) :: fmat(mat_dim, mat_dim, nz, nr_fmat)
      real(8), intent(in)    :: dmat_0(mat_dim, mat_dim, nz)
      real(8), intent(in)    :: dmat_0_saop_outer(mat_dim, mat_dim, nz)
      real(8), intent(in)    :: dmat_0_saop_resp(mat_dim, mat_dim, nz)
      integer, intent(in)    :: nr_dmat
      real(8), intent(in)    :: dmat_p(mat_dim, mat_dim, nz, nr_dmat)
      integer, intent(in)    :: dmat_pg_sym(nr_dmat)
      integer, intent(in)    :: dmat_ih_sym(nr_dmat)
      integer                :: derv_length
      real(8)                :: derv(max_block_length, 0:derv_length)
      integer, intent(in)    :: max_fun_derv
      real(8)                :: w(*)
      real(8)                :: nr_electrons_integrated
      real(8)                :: xc_energy
!     --------------------------------------------------------------------------
      real(8), allocatable   :: ao_compressed(:, :)
      integer                :: irep, imat
      integer                :: i, j, k, jcount
      integer                :: ixyz, jxyz
      integer                :: ikoff, jkoff, ijoff
      integer, parameter     :: max_nr_blocks = 16
      integer                :: isoff(max_nr_blocks)
      integer                :: iksoff(max_nr_blocks)
      integer                :: nr_above_threshold(max_nr_blocks)
      real(8)                :: temp
      integer, parameter     :: blas_limit = 10
      integer                :: is, js
      integer                :: nr1, nr2, st1, st2
      integer                :: iblock, jblock
      integer                :: jxyz_start
      real(8), allocatable   :: n(:, :)
      real(8), allocatable   :: n_saop_outer(:)
      real(8), allocatable   :: n_saop_resp(:)
      real(8), allocatable   :: n_g(:, :, :)
      real(8), allocatable   :: s(:, :)
      real(8), allocatable   :: v(:, :, :)
      real(8), allocatable   :: gnn(:, :, :)
      logical :: density_is_not_tiny 
!     --------------------------------------------------------------------------

      allocate(ao_compress_index(mat_dim))
      allocate(ao_compressed(max_block_length*mat_dim, 0:4))
      allocate(n(block_length, 0:nr_dmat))
      allocate(n_saop_outer(block_length))
      allocate(n_saop_resp(block_length))
      allocate(n_g(block_length, 3, 0:nr_dmat))
      allocate(s(block_length, maxval((/nr_dmat, 1/))))
      allocate(v(block_length, 3, maxval((/nr_dmat, 1/))))
      allocate(gnn(block_length, 0:1, 0:nr_dmat))

      density_is_not_tiny = .false.

!     zero out density and gradient
      do k = 1, block_length
         n(k, 0) = 0.0d0
      end do
      if (is_gga) then
         do ixyz = 1, 3
            do k = 1, block_length
               n_g(k, ixyz, 0) = 0.0d0
            end do
         end do
      end if
      if (dft_cfg_saop_is_active) then
         do k = 1, block_length
            n_saop_outer(k) = 0.0d0
         end do
         if (dft_cfg_saop_with_response_part) then
            do k = 1, block_length
               n_saop_resp(k) = 0.0d0
            end do
         end if
      end if

      isoff              = 0
      iksoff             = 0
      nr_above_threshold = 0

      if (nr_dmat > 0) then
         irep = dmat_pg_sym(1) - 1
      else
         irep = 0
      end if

      jcount = 0
      do iblock = 1, nr_ao_blocks
         st1 = ao_block_start(iblock)
         nr1 = ao_block_nr(iblock)

!        get buffer for compressing
!        and get number of functions above block_threshold
         do i = 1, nr1
            buffer(i) = 0.0d0
            is = i + st1 - 1
            do k = 1, block_length
               temp = dabs(ao_block(k, is))
               if (temp > buffer(i)) then
                  buffer(i) = temp
               end if
            end do
            if (buffer(i) > block_threshold) then
               jcount = jcount + 1
               ao_compress_index(jcount) = is
               nr_above_threshold(iblock) = nr_above_threshold(iblock) + 1
            end if
         end do

!        check for early return
         if (nr_above_threshold(iblock) == 0) then
            cycle
         end if

!        shift offsets of all later blocks
         do jblock = iblock + 1, max_nr_blocks
            isoff(jblock) = isoff(jblock) + nr_above_threshold(iblock)
         end do
         do jblock = iblock + 1, max_nr_blocks
            iksoff(jblock) = isoff(jblock)*block_length
         end do

!        compress ao matrix
         do i = 1, nr_above_threshold(iblock)
            is = ao_compress_index(isoff(iblock) + i)
            ikoff = (i - 1)*block_length
            do k = 1, block_length
               ao_compressed(iksoff(iblock) + ikoff + k, 0) = ao_block(k, is)
            end do
            if (is_gga) then
               do ixyz = 1, 3
                  do k = 1, block_length
                     ao_compressed(iksoff(iblock) + ikoff + k, ixyz) = ao_block(k, is + nr_ao_cartesian*ixyz)
                  end do
               end do
            end if
         end do
      end do

!     get unperturbed densities
      do iblock = 1, nr_ao_blocks
         if (nr_above_threshold(iblock) > 0) then
            call evaluate_density(                            &
                                  1.0d0,                      &
                                  is_gga,                     &
                                  block_length,               &
                                  mat_dim,                    &
                                  nz,                         &
                                  ao_compressed,                  &
                                  dmat_0,                     &
                                  nr_above_threshold(iblock), &
                                  nr_above_threshold(iblock), &
                                  isoff(iblock),              &
                                  isoff(iblock),              &
                                  iksoff(iblock),             &
                                  iksoff(iblock),             &
                                  n(1, 0),                    &
                                  n_g(1, 1, 0)                &
                                 )
            if (dft_cfg_saop_is_active) then
               call evaluate_density(                            &
                                     1.0d0,                      &
                                     .false.,                    &
                                     block_length,               &
                                     mat_dim,                    &
                                     nz,                         &
                                     ao_compressed,                  &
                                     dmat_0_saop_outer,          &
                                     nr_above_threshold(iblock), &
                                     nr_above_threshold(iblock), &
                                     isoff(iblock),              &
                                     isoff(iblock),              &
                                     iksoff(iblock),             &
                                     iksoff(iblock),             &
                                     n_saop_outer                &
                                    )
               if (dft_cfg_saop_with_response_part) then
                  call evaluate_density(                            &
                                        1.0d0,                      &
                                        .false.,                    &
                                        block_length,               &
                                        mat_dim,                    &
                                        nz,                         &
                                        ao_compressed,                  &
                                        dmat_0_saop_resp,           &
                                        nr_above_threshold(iblock), &
                                        nr_above_threshold(iblock), &
                                        isoff(iblock),              &
                                        isoff(iblock),              &
                                        iksoff(iblock),             &
                                        iksoff(iblock),             &
                                        n_saop_resp                 &
                                       )
               end if
            end if
         end if
      end do

      do k = 1, block_length
         if ((n(k,0).ge.dft_cfg_tinydens)) density_is_not_tiny = .true.
      end do

      if (dft_cfg_grac_is_active .or. dft_cfg_saop_is_active) then
         do k = 1, block_length
            if ((n_saop_resp(k).ge.dft_cfg_tinydens)) density_is_not_tiny = .true.
            if ((n_saop_outer(k).ge.dft_cfg_tinydens)) density_is_not_tiny = .true.
         end do
      end if

! if the density isn't below the small density cutoff we proceed to all
! evaluations below, otherwise we skip this block's contribution to the matrices

   if (density_is_not_tiny) then 
      if (is_gga) then
         do k = 1, block_length
            gnn(k, 0, 0) = n_g(k, 1, 0)*n_g(k, 1, 0) &
                         + n_g(k, 2, 0)*n_g(k, 2, 0) &
                         + n_g(k, 3, 0)*n_g(k, 3, 0)
         end do
      else
         do k = 1, block_length
            gnn(k, 0, 0) = 0.0d0
         end do
      end if

      call get_functional_derv(                           &
                               xc_fun,                    &
                               xc_fun_alda,               &
                               xc_fun_xalda,              &
                               max_fun_derv,              &
                               block_length,              &
                               w,                         &
                               n(1, 0),                   &
                               gnn(1, 0, 0),              &
                               derv_length,               &
                               derv,                      &
                               alda_real=dft_cfg_alda_hs, &
                               alda_imag=dft_cfg_alda_ha, &
                               xalda=dft_cfg_xalda        &
                              )

      if (dft_cfg_grac_is_active .or. dft_cfg_saop_is_active) then
         call dftac(block_length, derv, w, n, n_saop_resp, n_saop_outer, gnn)
      end if

      do k = 1, block_length
         nr_electrons_integrated = nr_electrons_integrated + w(k)*n(k, 0)
         xc_energy = xc_energy + derv(k, d0000000)
      end do

!     scf
      if (nr_dmat == 0) then
         do imat = 1, 1
!           form prefactors
            do k = 1, block_length
               s(k, imat) = derv(k, d1000000)
            end do
            if (is_gga) then
               do ixyz = 1, 3
                  do k = 1, block_length
                     v(k, ixyz, imat) = 2.0d0*n_g(k, ixyz, 0)*derv(k, d0010000)
                  end do
               end do
            end if
!           distribute
            do iblock = 1, nr_ao_blocks
               if (nr_above_threshold(iblock) > 0) then
                  call distribute_prefactors(                            &
                                             1.0d0,                      &
                                             is_gga,                     &
                                             block_length,               &
                                             mat_dim,                    &
                                             nz,                         &
                                             fmat(1, 1, 1, imat),        &
                                             ao_compressed,                  &
                                             nr_above_threshold(iblock), &
                                             isoff(iblock),              &
                                             iksoff(iblock),             &
                                             nr_above_threshold(iblock), &
                                             isoff(iblock),              &
                                             iksoff(iblock),             &
                                             s(1, imat),                 &
                                             v(1, 1, imat)               &
                                            )
               end if
            end do
         end do
      end if

      jxyz_start = 1
      if (dft_cfg_sdft_collinear) jxyz_start = 3 !collinear approximation

!     response
      if (nr_dmat > 0) then
         do imat = 1, nr_dmat
            irep = dmat_pg_sym(imat) - 1
            if (dmat_ih_sym(imat) == 1) then
               do jxyz = 0, 0
!                 get perturbed densities
                  do k = 1, block_length
                     n(k, imat) = 0.0d0
                  end do
                  if (is_gga_qr) then
                     do ixyz = 1, 3
                        do k = 1, block_length
                           n_g(k, ixyz, imat) = 0.0d0
                        end do
                     end do
                  end if
                  do iblock = 1, nr_ao_blocks
                     jblock = llss_block_partner(iblock, jxyz, irep)
                     if (jblock == 0) cycle
                     if ((nr_above_threshold(iblock) > 0) .and. (nr_above_threshold(jblock) > 0)) then
                        call evaluate_density(                              &
                                              1.0d0,                        &
                                              is_gga_qr,                    &
                                              block_length,                 &
                                              mat_dim,                      &
                                              nz,                           &
                                              ao_compressed,                    &
                                              dmat_p(1, 1, iq(jxyz), imat), &
                                              nr_above_threshold(iblock),   &
                                              nr_above_threshold(jblock),   &
                                              isoff(iblock),                &
                                              isoff(jblock),                &
                                              iksoff(iblock),               &
                                              iksoff(jblock),               &
                                              n(1, imat),                   &
                                              n_g(1, 1, imat)               &
                                             )
                     end if
                  end do
!                 form prefactors
                  do k = 1, block_length
                     s(k, imat) = derv(k, d2000000)*n(k, imat)
                  end do
                  if (is_gga_qr) then
                     do k = 1, block_length
                        gnn(k, 0, imat) = 2.0d0*n_g(k, 1, 0)*n_g(k, 1, imat) &
                                        + 2.0d0*n_g(k, 2, 0)*n_g(k, 2, imat) &
                                        + 2.0d0*n_g(k, 3, 0)*n_g(k, 3, imat)
                     end do
                     do k = 1, block_length
                        s(k, imat) = s(k, imat) + derv(k, d1010000)*gnn(k, 0, imat)
                     end do
                     do ixyz = 1, 3
                        do k = 1, block_length
                           v(k, ixyz, imat) = 2.0d0*n_g(k, ixyz, 0)*n(k, imat)*derv(k, d1010000)      &
                                            + 2.0d0*n_g(k, ixyz, 0)*gnn(k, 0, imat)*derv(k, d0020000) &
                                            + 2.0d0*n_g(k, ixyz, imat)*derv(k, d0010000)
                        end do
                     end do
                  end if
!                 distribute
                  do iblock = 1, nr_ao_blocks
                     jblock = llss_block_partner(iblock, jxyz, irep)
                     if (jblock == 0) cycle
                     if ((nr_above_threshold(iblock) > 0) .and. (nr_above_threshold(jblock) > 0)) then
                        call distribute_prefactors(                            &
                                                   1.0d0,                      &
                                                   is_gga_qr,                  &
                                                   block_length,               &
                                                   mat_dim,                    &
                                                   nz,                         &
                                                   fmat(1, 1, iq(jxyz), imat), &
                                                   ao_compressed,                  &
                                                   nr_above_threshold(iblock), &
                                                   isoff(iblock),              &
                                                   iksoff(iblock),             &
                                                   nr_above_threshold(jblock), &
                                                   isoff(jblock),              &
                                                   iksoff(jblock),             &
                                                   s(1, imat),                 &
                                                   v(1, 1, imat)               &
                                                  )
                     end if
                  end do
               end do
            end if
            if (dmat_ih_sym(imat) == -1 .and. .not. dft_cfg_no_sdft) then
               do jxyz = jxyz_start, 3
!                 get perturbed densities
                  do k = 1, block_length
                     n(k, imat) = 0.0d0
                  end do
                  if (is_gga_qi) then
                     do ixyz = 1, 3
                        do k = 1, block_length
                           n_g(k, ixyz, imat) = 0.0d0
                        end do
                     end do
                  end if
                  do iblock = 1, nr_ao_blocks
                     jblock = llss_block_partner(iblock, jxyz, irep)
                     if (jblock == 0) cycle
                     if ((nr_above_threshold(iblock) > 0) .and. (nr_above_threshold(jblock) > 0)) then
                        call evaluate_density(                              &
                                              llss_prefactor(iblock),       &
                                              is_gga_qi,                    &
                                              block_length,                 &
                                              mat_dim,                      &
                                              nz,                           &
                                              ao_compressed,                    &
                                              dmat_p(1, 1, iq(jxyz), imat), &
                                              nr_above_threshold(iblock),   &
                                              nr_above_threshold(jblock),   &
                                              isoff(iblock),                &
                                              isoff(jblock),                &
                                              iksoff(iblock),               &
                                              iksoff(jblock),               &
                                              n(1, imat),                   &
                                              n_g(1, 1, imat)               &
                                             )
                     end if
                  end do
!                 form prefactors
                  do k = 1, block_length
                     s(k, imat) = derv(k, d0200000)*n(k, imat)
                  end do
                  if (is_gga_qi) then
                     do k = 1, block_length
                        gnn(k, 0, imat) = n_g(k, 1, 0)*n_g(k, 1, imat) &
                                        + n_g(k, 2, 0)*n_g(k, 2, imat) &
                                        + n_g(k, 3, 0)*n_g(k, 3, imat)
                     end do
                     do k = 1, block_length
                        s(k, imat) = s(k, imat) + derv(k, d0101000)*gnn(k, 0, imat)
                     end do
                     do ixyz = 1, 3
                        do k = 1, block_length
                           v(k, ixyz, imat) =       n_g(k, ixyz, 0)*n(k, imat)*derv(k, d0101000)      &
                                            +       n_g(k, ixyz, 0)*gnn(k, 0, imat)*derv(k, d0002000) &
                                            + 2.0d0*n_g(k, ixyz, imat)*derv(k, d0000100)
                        end do
                     end do
                  end if
!                 distribute
                  do iblock = 1, nr_ao_blocks
                     jblock = llss_block_partner(iblock, jxyz, irep)
                     if (jblock == 0) cycle
                     if ((nr_above_threshold(iblock) > 0) .and. (nr_above_threshold(jblock) > 0)) then
                        call distribute_prefactors(                            &
                                                   llss_prefactor(iblock),     &
                                                   is_gga_qi,                  &
                                                   block_length,               &
                                                   mat_dim,                    &
                                                   nz,                         &
                                                   fmat(1, 1, iq(jxyz), imat), &
                                                   ao_compressed,                  &
                                                   nr_above_threshold(iblock), &
                                                   isoff(iblock),              &
                                                   iksoff(iblock),             &
                                                   nr_above_threshold(jblock), &
                                                   isoff(jblock),              &
                                                   iksoff(jblock),             &
                                                   s(1, imat),                 &
                                                   v(1, 1, imat)               &
                                                  )
                     end if
                  end do
               end do
            end if
         end do
      end if
  end if

      deallocate(ao_compress_index)
      deallocate(ao_compressed)
      deallocate(s)
      deallocate(v)
      deallocate(n)
      deallocate(n_saop_outer)
      deallocate(n_saop_resp)
      deallocate(n_g)
      deallocate(gnn)

   end subroutine

   subroutine distribute_prefactors(                      &
                                    f,                    &
                                    is_gga,               &
                                    block_length,         &
                                    mat_dim,              &
                                    nz,                   &
                                    fmat,                 &
                                    ao_compressed,            &
                                    nr_above_threshold_l, &
                                    isoff_l,              &
                                    iksoff_l,             &
                                    nr_above_threshold_r, &
                                    isoff_r,              &
                                    iksoff_r,             &
                                    s,                    &
                                    v                     &
                                   )

!     --------------------------------------------------------------------------
      real(8), intent(in)    :: f
      logical                :: is_gga
      integer                :: block_length
      integer                :: mat_dim
      integer                :: nz
      real(8)                :: ao_compressed(max_block_length*mat_dim, 0:4)
      real(8), intent(inout) :: fmat(mat_dim, mat_dim)
      integer                :: nr_above_threshold_l
      integer                :: isoff_l
      integer                :: iksoff_l
      integer                :: nr_above_threshold_r
      integer                :: isoff_r
      integer                :: iksoff_r
      real(8)                :: s(block_length)
      real(8)                :: v(block_length, 3)
!     --------------------------------------------------------------------------
      integer :: i, k, ixyz, j, is, ijoff, js, ikoff

      real(8), allocatable   :: xmat(:)
      real(8), allocatable   :: tmat(:)
!     --------------------------------------------------------------------------

      allocate(xmat(nr_above_threshold_l*nr_above_threshold_r))
      allocate(tmat(max_block_length*mat_dim))

      do i = 1, nr_above_threshold_l
         ikoff = (i - 1)*block_length
         do k = 1, block_length
            tmat(ikoff + k) = f*ao_compressed(iksoff_l + ikoff + k, 0)*s(k)
         end do
         if (is_gga) then
            do ixyz = 1, 3
               do k = 1, block_length
                  tmat(ikoff + k) = tmat(ikoff + k) + f*2.0d0*ao_compressed(iksoff_l + ikoff + k, ixyz)*v(k, ixyz)
               end do
            end do
         end if
      end do
      
      if (iksoff_l == iksoff_r) then !diagonal block
         call dsyr2k('u',                        & !uplo 
                     't',                        & !trans
                     nr_above_threshold_l,       & !n    
                     block_length,               & !k    
                     1.0d0,                      & !alpha
                     ao_compressed(iksoff_l + 1, 0), & !a    
                     block_length,               & !lda  
                     tmat,                       & !b    
                     block_length,               & !ldb  
                     0.0d0,                      & !beta 
                     xmat,                       & !c    
                     nr_above_threshold_l)         !ldc  
      else !off-diagonal block
         call dgemm('t',                         & !transa
                    'n',                         & !transb 
                    nr_above_threshold_r,        & !m      
                    nr_above_threshold_l,        & !n      
                    block_length,                & !k      
                    2.0d0,                       & !alpha  
                    ao_compressed(iksoff_r + 1, 0),  & !a    
                    block_length,                & !lda    
                    tmat,                        & !b      
                    block_length,                & !ldb    
                    0.0d0,                       & !beta   
                    xmat,                        & !c       
                    nr_above_threshold_r)          !ldc     
      end if
      
!     uncompress onto fmat
      if (iksoff_l == iksoff_r) then !diagonal block
         do i = 1, nr_above_threshold_l
            is = ao_compress_index(isoff_l + i)
            ijoff = (i - 1)*nr_above_threshold_l
!           diagonal
            fmat(is, is) = fmat(is, is) + 0.5d0*xmat(ijoff + i)
            do j = 1, i - 1
               js = ao_compress_index(isoff_l + j)
!              upper triangle
               fmat(js, is) = fmat(js, is) + xmat(ijoff + j)
            end do
         end do
      else !off-diagonal block
         do i = 1, nr_above_threshold_l
            is = ao_compress_index(isoff_l + i)
            ijoff = (i - 1)*nr_above_threshold_r
            do j = 1, nr_above_threshold_r
               js = ao_compress_index(isoff_r + j)
               fmat(js, is) = fmat(js, is) + 0.5d0*xmat(ijoff + j)
            end do
         end do
      end if

      deallocate(xmat)
      deallocate(tmat)

   end subroutine

   subroutine evaluate_density(                      &
                               f,                    &
                               is_gga,               &
                               block_length,         &
                               mat_dim,              &
                               nz,                   &
                               ao_compressed,            &
                               dmat,                 &
                               nr_above_threshold_l, &
                               nr_above_threshold_r, &
                               isoff_l,              &
                               isoff_r,              &
                               iksoff_l,             &
                               iksoff_r,             &
                               n,                    &
                               n_g                   &
                              )

!     --------------------------------------------------------------------------
      real(8), intent(in)              :: f
      logical, intent(in)              :: is_gga
      integer, intent(in)              :: block_length
      integer, intent(in)              :: mat_dim
      integer, intent(in)              :: nz
      real(8), intent(in)              :: ao_compressed(max_block_length*mat_dim, 0:4)
      real(8), intent(in)              :: dmat(mat_dim, mat_dim, nz)
      integer, intent(in)              :: nr_above_threshold_l
      integer, intent(in)              :: nr_above_threshold_r
      integer, intent(in)              :: isoff_l
      integer, intent(in)              :: isoff_r
      integer, intent(in)              :: iksoff_l
      integer, intent(in)              :: iksoff_r
      real(8), intent(inout)           :: n(block_length)
      real(8), intent(inout), optional :: n_g(block_length, 3)
!     --------------------------------------------------------------------------
      integer                          :: ixyz
      real(8), allocatable             :: dmat_compressed(:)
      real(8), allocatable             :: xmat(:)
!     --------------------------------------------------------------------------

      allocate(xmat(mat_dim*max_block_length))

      allocate(dmat_compressed(nr_above_threshold_l*nr_above_threshold_r))

      call compress_mat(                      &
                    mat_dim,              &
                    dmat,                 &
                    dmat_compressed,          &
                    nr_above_threshold_l, &
                    nr_above_threshold_r, &
                    isoff_l,              &
                    isoff_r               &
                   )

      if (iksoff_l == iksoff_r) then
         call dsymm(                            &
                    'r',                        & !side   
                    'u',                        & !uplo   
                    block_length,               & !m      
                    nr_above_threshold_l,       & !n      
                    1.0d0,                      & !alpha  
                    dmat_compressed,                & !a      
                    nr_above_threshold_l,       & !lda    
                    ao_compressed(iksoff_l + 1, 0), & !b      
                    block_length,               & !ldb    
                    0.0d0,                      & !beta   
                    xmat,                       & !c      
                    block_length                & !ldc     
                   )
      else
         call dgemm(                            &
                    'n',                        & !transa
                    'n',                        & !transb 
                    block_length,               & !m      
                    nr_above_threshold_r,       & !n      
                    nr_above_threshold_l,       & !k      
                    1.0d0,                      & !alpha  
                    ao_compressed(iksoff_l + 1, 0), & !a      
                    block_length,               & !lda    
                    dmat_compressed,                & !b    
                    nr_above_threshold_l,       & !ldb
                    0.0d0,                      & !beta   
                    xmat,                       & !c       
                    block_length                & !ldc     
                   )
      end if

      call multiply_xmat_with_ao(                      &
                                 f,                    &
                                 block_length,         &
                                 nr_above_threshold_r, &
                                 iksoff_r,             &
                                 xmat,                 &
                                 ao_compressed(1, 0),      &
                                 n                     &
                                )
      if (is_gga) then
         do ixyz = 1, 3
            call multiply_xmat_with_ao(                      &
                                       2.0d0*f,              &
                                       block_length,         &
                                       nr_above_threshold_r, &
                                       iksoff_r,             &
                                       xmat,                 &
                                       ao_compressed(1, ixyz),   &
                                       n_g(1, ixyz)          &
                                      )
         end do
      end if

      deallocate(dmat_compressed)
      deallocate(xmat)

   end subroutine

   function get_block_threshold(         &
                                dfthr0,  &
                                mat_dim, &
                                dmat_0   &
                               )


!     --------------------------------------------------------------------------
      real(8)             :: get_block_threshold
!     --------------------------------------------------------------------------
      real(8), intent(in) :: dfthr0
      integer, intent(in) :: mat_dim
      real(8), intent(in) :: dmat_0(mat_dim, mat_dim)
!     --------------------------------------------------------------------------
      integer             :: i, j, iblock, st, nr
      real(8)             :: f, x
!     --------------------------------------------------------------------------

!     x = 0.0d0
!     iblock_loop: do iblock = 1, nr_ao_blocks
!        st = ao_block_start(iblock)
!        nr = ao_block_nr(iblock)
!        do i = st, st + nr - 1
!           do j = st, i
!              dmat_0 is positive definite but just to be sure
!              f = dabs(dmat_0(j, i))
!              if (f > x) then
!                 x = f
!              end if
!           end do
!        end do
!     end do iblock_loop
!radovan: my own empirical estimate, no physics behind this
!     get_block_threshold = dsqrt(10.0d0*dfthr0/(mat_dim*mat_dim*x))
      get_block_threshold = dft_cfg_screening

   end function

   subroutine compress_mat(                      &
                       mat_dim,              &
                       mat,                  &
                       mat_compressed,           &
                       nr_above_threshold_l, &
                       nr_above_threshold_r, &
                       isoff_l,              &
                       isoff_r               &
                      )

!     --------------------------------------------------------------------------
      integer, intent(in)  :: mat_dim
      real(8), intent(in)  :: mat(*)
      real(8), intent(out) :: mat_compressed(*)
      integer, intent(in)  :: nr_above_threshold_l
      integer, intent(in)  :: nr_above_threshold_r
      integer, intent(in)  :: isoff_l
      integer, intent(in)  :: isoff_r
!     --------------------------------------------------------------------------
      integer              :: ijoff
      integer              :: i, j
      integer              :: is, js
!     --------------------------------------------------------------------------

      do i = 1, nr_above_threshold_r
         is = ao_compress_index(isoff_r + i)
         ijoff = (i - 1)*nr_above_threshold_l
         do j = 1, nr_above_threshold_l
            js = ao_compress_index(isoff_l + j)
            mat_compressed(ijoff + j) = mat(mat_dim*(is - 1) + js)
         end do
      end do

   end subroutine

   subroutine multiply_xmat_with_ao(                    &
                                    f,                  &
                                    block_length,       &
                                    nr_above_threshold, &
                                    iksoff,             &
                                    xmat,               &
                                    ao_compressed,          &
                                    density             &
                                   )

!     --------------------------------------------------------------------------
      real(8), intent(in)    :: f
      integer, intent(in)    :: block_length
      integer, intent(in)    :: nr_above_threshold
      integer, intent(in)    :: iksoff
      real(8), intent(in)    :: xmat(*)
      real(8), intent(in)    :: ao_compressed(*)
      real(8), intent(inout) :: density(*)
!     --------------------------------------------------------------------------
      integer                :: i, k, ikoff
!     --------------------------------------------------------------------------

      do i = 1, nr_above_threshold
         ikoff = (i - 1)*block_length
         do k = 1, block_length
            density(k) = density(k) + f*xmat(ikoff + k)*ao_compressed(iksoff + ikoff + k)
         end do
      end do

   end subroutine

end module
