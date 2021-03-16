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

module xc_london_c1

! radovan: this is the "equation on paper approach" to LAO XC
!
!      (+) the code is very close to equations
!      (+) code does not need "magnetic derivatives" of AOs
!          and works with undifferentiated pointwise AO overlap
!      (+) i think it should be very easy to do higher order LAO
!          contributions to KS matrix and to KS response
!          using this approach
!          it should boil down to additional factorization of what we already have
!      (-) it does not work in symmetry adapted (SA) AO basis
!          because we cannot calculate SA AOs and "SA" R_kl in an uncoupled way
!          one could still use R_kl but introduce additional loops
!          in distribution routines and work with non-SA AOs
!          and transform to SA in each matrix element
!          but i think this would be very expensive
!          anyway this code can be useful for development in C1

! todo (*) work only on matrix blocks that couple AOs of different centers

   use interface_ao
   use dirac_ao_eval

   implicit none

   public xc_london_c1_ks_lda
   public xc_london_c1_ks_gga
   public xc_london_c1_lr_lda_densities
   public xc_london_c1_lr_gga_densities
   public xc_london_c1_get_s_b
   public xc_london_c1_get_gs_b
   public xc_london_c1_susc_der1
   public xc_london_c1_susc_der2
   public j_dia
   public xc_london_c1_susc_gga_der1

   private

   integer :: b(3) = (/2, 3, 1/)
   integer :: c(3) = (/3, 1, 2/)

contains

   subroutine xc_london_c1_ks_lda(d_n,     &
                                  mat_dim, &
                                  nz,      &
                                  fmat,    &
                                  ao,      &
                                  r)

!     --------------------------------------------------------------------------
      real(8), intent(in)    :: d_n
      integer, intent(in)    :: mat_dim
      integer, intent(in)    :: nz
      real(8), intent(inout) :: fmat(mat_dim, mat_dim, *)
      real(8), intent(in)    :: ao(mat_dim)
      real(8), intent(in)    :: r(3)
!     --------------------------------------------------------------------------
      integer                :: k, l, bxyz
      real(8)                :: f
      real(8)                :: rkl_cross_r
      real(8)                :: omega
!     --------------------------------------------------------------------------

      do k = 1, mat_dim
!        work only on the triangle without diagonal elements
         do l = 1, k - 1

            omega = ao(k)*ao(l)

            do bxyz = 1, 3
               rkl_cross_r = ao_distance(k, l, b(bxyz))*r(c(bxyz)) &
                           - ao_distance(k, l, c(bxyz))*r(b(bxyz))
               rkl_cross_r = 0.5d0*rkl_cross_r

               f = d_n*rkl_cross_r*omega

               fmat(k, l, (bxyz-1)*nz + 1) = fmat(k, l, (bxyz-1)*nz + 1) + f
               fmat(l, k, (bxyz-1)*nz + 1) = fmat(l, k, (bxyz-1)*nz + 1) - f
            end do
         end do
      end do

   end subroutine

   subroutine xc_london_c1_ks_gga(d_n,     &
                                  d_gnn,   &
                                  gn_0,    &
                                  mat_dim, &
                                  nz,      &
                                  fmat,    &
                                  ao,      &
                                  r)

!     --------------------------------------------------------------------------
      real(8), intent(in)    :: d_n
      real(8), intent(in)    :: d_gnn
      real(8), intent(in)    :: gn_0(3)
      integer, intent(in)    :: mat_dim
      integer, intent(in)    :: nz
      real(8), intent(inout) :: fmat(mat_dim, mat_dim, *)
      real(8), intent(in)    :: ao(*)
      real(8), intent(in)    :: r(3)
!     --------------------------------------------------------------------------
      integer                :: k, l, bxyz, gxyz
      real(8)                :: f, gnn_kl
      real(8)                :: rkl_cross_r, rkl_cross_gn_0
      real(8)                :: omega, nabla_omega(3)
!     --------------------------------------------------------------------------

      do k = 1, mat_dim
!        work only on the triangle without diagonal elements
         do l = 1, k - 1

            omega = ao(k)*ao(l)
            do gxyz = 1, 3
               nabla_omega(gxyz) = ao(k)*ao(ao_off_g1_m0(gxyz, 0) + l) &
                                 + ao(l)*ao(ao_off_g1_m0(gxyz, 0) + k)
            end do

            gnn_kl = gn_0(1)*nabla_omega(1) &
                   + gn_0(2)*nabla_omega(2) &
                   + gn_0(3)*nabla_omega(3)
            gnn_kl = 2.0d0*gnn_kl

            do bxyz = 1, 3
               rkl_cross_r = ao_distance(k, l, b(bxyz))*r(c(bxyz)) &
                           - ao_distance(k, l, c(bxyz))*r(b(bxyz))
               rkl_cross_r = 0.5d0*rkl_cross_r

               rkl_cross_gn_0 = ao_distance(k, l, b(bxyz))*gn_0(c(bxyz)) &
                              - ao_distance(k, l, c(bxyz))*gn_0(b(bxyz))
               rkl_cross_gn_0 = 0.5d0*rkl_cross_gn_0

               f = d_n*rkl_cross_r*omega    &
                 + d_gnn*rkl_cross_r*gnn_kl &
                 + d_gnn*2.0d0*rkl_cross_gn_0*omega

               fmat(k, l, (bxyz-1)*nz + 1) = fmat(k, l, (bxyz-1)*nz + 1) + f
               fmat(l, k, (bxyz-1)*nz + 1) = fmat(l, k, (bxyz-1)*nz + 1) - f
            end do
         end do
      end do

   end subroutine

   subroutine xc_london_c1_lr_lda_densities(mat_dim, &
                                            nz,      &
                                            dmat,    &
                                            ao,      &
                                            r,       &
                                            n_b)

!     --------------------------------------------------------------------------
      integer, intent(in)  :: mat_dim
      integer, intent(in)  :: nz
      real(8), intent(in)  :: dmat(mat_dim, mat_dim)
      real(8), intent(in)  :: ao(mat_dim)
      real(8), intent(in)  :: r(3)
      real(8), intent(out) :: n_b(3)
!     --------------------------------------------------------------------------
      integer              :: k, l, bxyz
      real(8)              :: f
      real(8)              :: rkl_cross_r
      real(8)              :: omega
!     --------------------------------------------------------------------------

      n_b = 0.0d0

      do k = 1, mat_dim
         do l = 1, mat_dim

            f = -1.0d0*dmat(l, k)
            omega = ao(k)*ao(l)*f

            do bxyz = 1, 3
               rkl_cross_r = ao_distance(k, l, b(bxyz))*r(c(bxyz)) &
                           - ao_distance(k, l, c(bxyz))*r(b(bxyz))
               rkl_cross_r = 0.5d0*rkl_cross_r

               n_b(bxyz) = n_b(bxyz) &
                         + rkl_cross_r*omega
            end do
         end do
      end do

   end subroutine

   subroutine xc_london_c1_lr_gga_densities(mat_dim, &
                                            nz,      &
                                            dmat,    &
                                            ao,      &
                                            r,       &
                                            gn_0,    &
                                            n_b,     &
                                            gn_b,    &
                                            gnn_b)

!     --------------------------------------------------------------------------
      integer, intent(in)  :: mat_dim
      integer, intent(in)  :: nz
      real(8), intent(in)  :: dmat(mat_dim, mat_dim)
      real(8), intent(in)  :: ao(*)
      real(8), intent(in)  :: r(3)
      real(8), intent(in)  :: gn_0(3)
      real(8), intent(out) :: n_b(3)
      real(8), intent(out) :: gn_b(3, 3)
      real(8), intent(out) :: gnn_b(3)
!     --------------------------------------------------------------------------
      integer              :: k, l, bxyz, gxyz
      real(8)              :: f, gnn_kl
      real(8)              :: rkl_cross_r, rkl_cross_gn_0
      real(8)              :: omega, nabla_omega(3)
!     --------------------------------------------------------------------------

      n_b   = 0.0d0
      gn_b  = 0.0d0
      gnn_b = 0.0d0

      do k = 1, mat_dim
         do l = 1, mat_dim

            f = -1.0d0*dmat(l, k)
            omega = ao(k)*ao(l)*f
            do gxyz = 1, 3
               nabla_omega(gxyz) = ao(k)*ao(ao_off_g1_m0(gxyz, 0) + l) &
                                 + ao(l)*ao(ao_off_g1_m0(gxyz, 0) + k)
               nabla_omega(gxyz) = f*nabla_omega(gxyz)
            end do

            gnn_kl = gn_0(1)*nabla_omega(1) &
                   + gn_0(2)*nabla_omega(2) &
                   + gn_0(3)*nabla_omega(3)
            gnn_kl = 2.0d0*gnn_kl

            do bxyz = 1, 3
               rkl_cross_r = ao_distance(k, l, b(bxyz))*r(c(bxyz)) &
                           - ao_distance(k, l, c(bxyz))*r(b(bxyz))
               rkl_cross_r = 0.5d0*rkl_cross_r

               rkl_cross_gn_0 = ao_distance(k, l, b(bxyz))*gn_0(c(bxyz)) &
                              - ao_distance(k, l, c(bxyz))*gn_0(b(bxyz))
               rkl_cross_gn_0 = 0.5d0*rkl_cross_gn_0

               n_b(bxyz) = n_b(bxyz) &
                         + rkl_cross_r*omega

               gnn_b(bxyz) = gnn_b(bxyz)        &
                           + rkl_cross_r*gnn_kl &
                           + 2.0d0*rkl_cross_gn_0*omega

               do gxyz = 1, 3
                  gn_b(gxyz, bxyz) = gn_b(gxyz, bxyz)              &
                                   + rkl_cross_r*nabla_omega(gxyz)
               end do
            end do
!           gn_b(1, 1)
            gn_b(1, 2) = gn_b(1, 2) + 0.5d0*ao_distance(k, l, 3)*omega
            gn_b(1, 3) = gn_b(1, 3) - 0.5d0*ao_distance(k, l, 2)*omega
            gn_b(2, 1) = gn_b(2, 1) - 0.5d0*ao_distance(k, l, 3)*omega
!           gn_b(2, 2)
            gn_b(2, 3) = gn_b(2, 3) + 0.5d0*ao_distance(k, l, 1)*omega
            gn_b(3, 1) = gn_b(3, 1) + 0.5d0*ao_distance(k, l, 2)*omega
            gn_b(3, 2) = gn_b(3, 2) - 0.5d0*ao_distance(k, l, 1)*omega
!           gn_b(3, 3)
         end do
      end do

   end subroutine

   subroutine j_dia(j,            &
                    lao,          &
                    gauge_origin, &
                    irep,         &
                    mat_dim,      &
                    dmat,         &
                    ao,           &
                    r)

!     --------------------------------------------------------------------------
      real(8), intent(out)   :: j(3, 3)
      logical, intent(in)    :: lao
      real(8), intent(in)    :: gauge_origin(3)
      integer, intent(in)    :: irep
      integer, intent(in)    :: mat_dim
      real(8), intent(inout) :: dmat(mat_dim, mat_dim, *)
      real(8), intent(in)    :: ao(mat_dim)
      real(8), intent(in)    :: r(3)
!     --------------------------------------------------------------------------
      integer                :: k, l, ixyz
      real(8)                :: f
      real(8)                :: q(3)
      integer                :: iblock, nr1, nr2, st1, st2
!     --------------------------------------------------------------------------

      if (lao .and. nr_boson_ireps > 1) then
         print *, 'ERROR: lao diamagnetic current only works in C1 sym'
         stop 1
      end if

      q = 0.0d0

      do iblock = 1, nr_ao_blocks
         nr1 = ao_block_nr(iblock)
         st1 = ao_block_start(iblock)
         nr2 = ao_block_nr(llss_block_partner(iblock, 0, irep))
         st2 = ao_block_start(llss_block_partner(iblock, 0, irep))
         if ((nr1 < 1) .or. (nr2 < 1)) cycle

         if (lao) then
            do k = st1, st1 + nr1 - 1
               do l = st2, st2 + nr2 - 1
                  f = ao(k)*ao(l)*dmat(l, k, 1)
                  do ixyz = 1, 3
                     q(ixyz) = q(ixyz) + (r(ixyz) - ao_center_xyz(l, ixyz))*f
                  end do
               end do
            end do
         else
            do k = st1, st1 + nr1 - 1
               do l = st2, st2 + nr2 - 1
                  f = ao(k)*ao(l)*dmat(l, k, 1)
                  do ixyz = 1, 3
                     q(ixyz) = q(ixyz) + (r(ixyz) - gauge_origin(ixyz))*f
                  end do
               end do
            end do
         end if
      end do

      j = 0.0d0

      j(1, 2) =  q(3) ! d j_x / d B_y
      j(2, 1) = -q(3) ! d j_y / d B_x

      j(3, 1) =  q(2) ! d j_z / d B_x
      j(1, 3) = -q(2) ! d j_x / d B_z

      j(2, 3) =  q(1) ! d j_y / d B_z
      j(3, 2) = -q(1) ! d j_z / d B_y

      j = 0.5d0*j

   end subroutine

   subroutine xc_london_c1_get_s_b(s,       &
                                   mat_dim, &
                                   nz,      &
                                   dmat,    &
                                   ao,      &
                                   r)

!     --------------------------------------------------------------------------
      real(8), intent(out) :: s(3, 3)
      integer, intent(in)  :: mat_dim
      integer, intent(in)  :: nz
      real(8), intent(in)  :: dmat(mat_dim, mat_dim, nz)
      real(8), intent(in)  :: ao(mat_dim)
      real(8), intent(in)  :: r(3)
!     --------------------------------------------------------------------------
      integer              :: k, l, n, bxyz
      real(8)              :: f
      real(8)              :: rkl_cross_r
      real(8)              :: omega
!     --------------------------------------------------------------------------

      s = 0.0d0

      do n = 1, 3
         do k = 1, mat_dim
            do l = 1, mat_dim
         
               f = dmat(l, k, iq(n))
               omega = ao(k)*ao(l)*f
         
               do bxyz = 1, 3
                  rkl_cross_r = ao_distance(k, l, b(bxyz))*r(c(bxyz)) &
                              - ao_distance(k, l, c(bxyz))*r(b(bxyz))
                  rkl_cross_r = 0.5d0*rkl_cross_r
         
                  s(n, bxyz) = s(n, bxyz) &
                            + rkl_cross_r*omega
               end do
            end do
         end do
      end do

   end subroutine

   subroutine xc_london_c1_get_gs_b(gs,      &
                                    s,       &
                                    mat_dim, &
                                    nz,      &
                                    dmat,    &
                                    ao,      &
                                    r)

!     --------------------------------------------------------------------------
      real(8), intent(out) :: gs(3, 3, 3)
      real(8), intent(out) :: s(3, 3)
      integer, intent(in)  :: mat_dim
      integer, intent(in)  :: nz
      real(8), intent(in)  :: dmat(mat_dim, mat_dim, nz)
      real(8), intent(in)  :: ao(mat_dim)
      real(8), intent(in)  :: r(3)
!     --------------------------------------------------------------------------
      integer              :: k, l, n, bxyz, ixyz
      real(8)              :: f
      real(8)              :: rkl_cross_r, nabla_rkl_cross_r
      real(8)              :: rkl_cross_gn
      real(8)              :: gn_nabla_omega
      real(8)              :: omega, nabla_omega(3)
!     --------------------------------------------------------------------------

      s = 0.0d0   ! \Sigma, B
      gs = 0.0d0  ! grad, \sigma, B

      do n = 1, 3   ! \Sigma_{x,y,z}
         do k = 1, mat_dim
            do l = 1, mat_dim
         
               f = dmat(l, k, iq(n))
               omega = ao(k)*ao(l)*f
         
               do ixyz = 1, 3 ! dR_{x,y,z}
                  nabla_omega(ixyz) = ao(k)*ao(ao_off_g1_m0(ixyz, 0) + l) &
                                    + ao(l)*ao(ao_off_g1_m0(ixyz, 0) + k)
                  nabla_omega(ixyz) = nabla_omega(ixyz)*f
               end do

               do bxyz = 1, 3  ! dB_{x,y,z}
                  rkl_cross_r = ao_distance(k, l, b(bxyz))*r(c(bxyz)) &
                              - ao_distance(k, l, c(bxyz))*r(b(bxyz))
                  rkl_cross_r = 0.5d0*rkl_cross_r
         
                  s(n, bxyz) = s(n, bxyz) &
                            + rkl_cross_r*omega

                  do ixyz = 1, 3

                     if (ixyz .eq. b(bxyz)) then
                         nabla_rkl_cross_r = - 0.5d0*ao_distance(k, l, c(bxyz))
                     else if (ixyz .eq. c(bxyz)) then
                         nabla_rkl_cross_r =   0.5d0*ao_distance(k, l, b(bxyz))
                     else
                        nabla_rkl_cross_r = 0.0d0
                     end if

                     gs(ixyz, n, bxyz) = gs(ixyz, n, bxyz) &
                                       + rkl_cross_r*nabla_omega(ixyz) &
                                       + nabla_rkl_cross_r*omega

                  end do
               end do

            end do
         end do
      end do

   end subroutine


   subroutine xc_london_c1_susc_der1(s, &
                                     mat_dim, &
                                     nz,      &
                                     dmat,    &
                                     ao,      &
                                     r,       &
                                     n_b)

!     --------------------------------------------------------------------------
      integer, intent(in)  :: mat_dim
      integer, intent(in)  :: nz
      real(8), intent(in)  :: s
      real(8), intent(in)  :: dmat(mat_dim, mat_dim)
      real(8), intent(in)  :: ao(mat_dim)
      real(8), intent(in)  :: r(3)
      real(8), intent(out) :: n_b(3, 3)
!     --------------------------------------------------------------------------
      integer              :: k, l, alpha, beta
      real(8)              :: f, f1, f2
      real(8)              :: rkl_cross_r(2)
      real(8)              :: omega
!     --------------------------------------------------------------------------

      n_b = 0.0d0

      do k = 1, mat_dim
         do l = 1, mat_dim

            f = -0.25d0*dmat(l, k)
            omega = ao(k)*ao(l)*f*s

            do alpha = 1, 3
              do beta = 1, 3

               f1 = ao_distance(k, l, b(alpha))*r(c(alpha)) &
                                  - ao_distance(k, l, c(alpha))*r(b(alpha))

               f2 = ao_distance(k, l, b(beta))*r(c(beta)) &
                                 - ao_distance(k, l, c(beta))*r(b(beta))

               n_b(beta, alpha) = n_b(beta, alpha) &
                         + f1*f2*omega

              end do
            end do
         end do
      end do


   end subroutine


   subroutine xc_london_c1_susc_gga_der1(d_r, &
                                     d_z,     &
                                     gn_0,    &
                                     mat_dim, &
                                     nz,      &
                                     dmat,    &
                                     ao,      &
                                     r,       &
                                     n_b)

!     --------------------------------------------------------------------------
      integer, intent(in)  :: mat_dim
      integer, intent(in)  :: nz
      real(8), intent(in)  :: d_r, d_z, gn_0(3)
      real(8), intent(in)  :: dmat(mat_dim, mat_dim)
      real(8), intent(in)  :: ao(mat_dim)
      real(8), intent(in)  :: r(3)
      real(8), intent(out) :: n_b(3, 3)
!     --------------------------------------------------------------------------
      integer              :: k, l, alpha, beta, gxyz
      real(8)              :: f, f1, f2, f3, f4, gnn_kl
      real(8)              :: omega, nabla_omega(3)
!     --------------------------------------------------------------------------

      n_b = 0.0d0
      gnn_kl = 0.0d0

      do k = 1, mat_dim
         do l = 1, mat_dim

            f = -0.25d0*dmat(l, k)
            omega = ao(k)*ao(l)*f

            do gxyz = 1, 3
               nabla_omega(gxyz) = ao(k)*ao(ao_off_g1_m0(gxyz, 0) + l) &
                                 + ao(l)*ao(ao_off_g1_m0(gxyz, 0) + k)
               nabla_omega(gxyz) = f*nabla_omega(gxyz)
            end do

            gnn_kl = gn_0(1)*nabla_omega(1) &
                   + gn_0(2)*nabla_omega(2) &
                   + gn_0(3)*nabla_omega(3)

            do alpha = 1, 3
              do beta = 1, 3

               f1 = ao_distance(k, l, b(alpha))*r(c(alpha)) &
                  - ao_distance(k, l, c(alpha))*r(b(alpha))

               f2 = ao_distance(k, l, b(beta))*r(c(beta)) &
                  - ao_distance(k, l, c(beta))*r(b(beta))

               f3 = ao_distance(k, l, b(alpha))*gn_0(c(alpha)) &
                  - ao_distance(k, l, c(alpha))*gn_0(b(alpha))

               f4 = ao_distance(k, l, b(beta))*gn_0(c(beta)) &
                  - ao_distance(k, l, c(beta))*gn_0(b(beta))


               n_b(beta, alpha) = n_b(beta, alpha) &
                         + f1*f2*omega*d_r  + (f1*f2*gnn_kl + f3*f2*omega + f4*f1*omega)*d_z

              end do
            end do
         end do
      end do

   end subroutine

   subroutine xc_london_c1_susc_der2(s, &
                                     mat_dim, &
                                     nz,      &
                                     dmat,    &
                                     ao,      &
                                     r,       &
                                     n_b)

!     --------------------------------------------------------------------------
      integer, intent(in)  :: mat_dim
      integer, intent(in)  :: nz
      real(8), intent(in)  :: s
      real(8), intent(in)  :: dmat(mat_dim, mat_dim)
      real(8), intent(in)  :: ao(mat_dim)
      real(8), intent(in)  :: r(3)
      real(8), intent(out) :: n_b(3, 3)
!     --------------------------------------------------------------------------
      integer              :: k, l, alpha, beta
      real(8)              :: f, f1, n1(3)
      real(8)              :: omega
!     --------------------------------------------------------------------------

      n_b = 0.0d0

      do k = 1, mat_dim
         do l = 1, mat_dim

            f = dmat(l, k)
            omega = ao(k)*ao(l)*f

            do alpha = 1, 3
               f1 = ao_distance(k, l, b(alpha))*r(c(alpha)) &
                           - ao_distance(k, l, c(alpha))*r(b(alpha))

               n1(alpha) = n1(alpha) + f1*omega
            end do
         end do
      end do

      do k = 1, 3
        do l = 1, 3
          n_b(l, k) = n_b(l, k) - 0.25d0*n1(l)*n1(k)*s
        end do
      end do


   end subroutine

end module
