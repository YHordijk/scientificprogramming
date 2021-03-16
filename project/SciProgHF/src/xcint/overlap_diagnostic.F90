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

module overlap_diagnostic

!  MO overlap diagnostic according to Peach, Benfield, Helgaker, and Tozer:
!  M. J. G. Peach, P. Benfield, T. U. Helgaker, and D. J. Tozer
!  J. Chem. Phys. 128, 044118 (2008)

!  implemented 2011/01/11 by Radovan Bast <radovan.bast@uit.no>

   use interface_ao
   use interface_mo
   use xc_mpi
   use xc_max_block_length

   implicit none

   public sum_lambda
   public integrate_o_matrix
   public overlap_diagnostic_init
   public overlap_diagnostic_collect
   public overlap_diagnostic_delete
   public o_matrix_to_disc
   public o_matrix_from_disc

   private

   integer                 :: nr_g_ao_blocks(2)
   integer                 :: nr_u_ao_blocks(2)
   integer                 :: g_ao_start(8,2)
   integer                 :: u_ao_start(8,2)
   integer                 :: g_ao_end(8,2)
   integer                 :: u_ao_end(8,2)
   integer                 :: nc
   real(8), allocatable,public    :: o_matrix(:,:)
   real(8), allocatable    :: spinor(:, :, :, :)

   character(*), parameter :: file_name = 'o_matrix'
   integer,      parameter :: file_unit = 137
   logical                 :: file_exists

contains

   subroutine sum_lambda(lambda_numerator,   &
                         lambda_denominator, &
                         i,                  &
                         s,                  &
                         ifermi,             &
                         sfermi,             &
                         kappa)

!     --------------------------------------------------------------------------
      real(8), intent(inout) :: lambda_numerator
      real(8), intent(inout) :: lambda_denominator
      integer, intent(in)    :: i
      integer, intent(in)    :: s
      integer, intent(in)    :: ifermi
      integer, intent(in)    :: sfermi
      real(8), intent(in)    :: kappa
!     --------------------------------------------------------------------------
      integer                :: ioff, soff
!     --------------------------------------------------------------------------

      select case (ifermi)
         case (1)
            ioff = nr_mo_gerade_negative_secondary
         case (2)
            ioff = nr_mo_gerade + nr_mo_ungerade_negative_secondary
      end select

      select case (sfermi)
         case (1)
            soff = nr_mo_gerade_negative_secondary
         case (2)
            soff = nr_mo_gerade + nr_mo_ungerade_negative_secondary
      end select

      lambda_numerator   = lambda_numerator   + kappa*o_matrix(ioff + i, soff + s)
      lambda_denominator = lambda_denominator + kappa

   end subroutine

   subroutine overlap_diagnostic_init(parallel_xc, i_am_master)


!     --------------------------------------------------------------------------
      logical, intent(in) :: parallel_xc, i_am_master
!     --------------------------------------------------------------------------
      integer             :: k, irep, ic
!     --------------------------------------------------------------------------

      if (i_am_master) then
         
!         nr_g_ao_blocks(1) = 0
!         nr_u_ao_blocks(1) = 0

         nr_g_ao_blocks = 0
         nr_u_ao_blocks = 0
         g_ao_start     = 0
         u_ao_start     = 0
         g_ao_end       = 0
         u_ao_end       = 0
         
!        fixme doing too much for Levy-Leblond and .X2C4
         if (nr_shells_small > 0) then
            nc = 2
         else
            nc = 1
           nr_g_ao_blocks(2) = 0
           nr_u_ao_blocks(2) = 0
         end if

         do ic = 1, nc
            do irep = 0, nr_boson_ireps - 1
               if (ao_nr(irep, ic) > 0) then
!                 fixme check this carefully
                  if (ao_fermion_corep(irep, ic) == 1) then
                     nr_g_ao_blocks(ic) = nr_g_ao_blocks(ic) + 1
!                 separating small and large component            
                     g_ao_start(nr_g_ao_blocks(ic),ic) = ao_off(irep, ic) + 1
                     g_ao_end(nr_g_ao_blocks(ic),ic)   = ao_off(irep, ic) + ao_nr(irep, ic)
                  else
                     nr_u_ao_blocks(ic) = nr_u_ao_blocks(ic) + 1
                     u_ao_start(nr_u_ao_blocks(ic),ic) = ao_off(irep, ic) + 1
                     u_ao_end(nr_u_ao_blocks(ic),ic)   = ao_off(irep, ic) + ao_nr(irep, ic)
                  end if
               end if
            end do
         end do
      end if

#ifdef VAR_MPI
      if (parallel_xc) then
         call xc_mpi_bcast(nc)
         call xc_mpi_bcast(nr_g_ao_blocks)
         call xc_mpi_bcast(nr_u_ao_blocks)
         call xc_mpi_bcast(g_ao_start    )
         call xc_mpi_bcast(u_ao_start    )
         call xc_mpi_bcast(g_ao_end      )
         call xc_mpi_bcast(u_ao_end      )
      end if
#endif

      allocate(o_matrix(nr_mo, nr_mo))
      allocate(spinor(max_block_length, nr_mo, 4,nc))

      o_matrix = 0.0d0

   end subroutine

   subroutine overlap_diagnostic_collect(i_am_master)


!     --------------------------------------------------------------------------
      logical             :: i_am_master
!     --------------------------------------------------------------------------
      integer             :: ierr
      real(8)             :: dummy
!     --------------------------------------------------------------------------

#ifdef VAR_MPI
         call xc_mpi_reduce(o_matrix)
#endif

   end subroutine

   subroutine overlap_diagnostic_delete()

      nr_g_ao_blocks = 0
      nr_u_ao_blocks = 0
      g_ao_start     = 0
      u_ao_start     = 0
      g_ao_end       = 0
      u_ao_end       = 0

      if (allocated(o_matrix)) deallocate(o_matrix)
      if (allocated(spinor))   deallocate(spinor)

   end subroutine

   subroutine o_matrix_to_disc()

      inquire(file = file_name, exist = file_exists)

      if (file_exists) then
!        panic
!        close(file_unit, status = 'keep')
!        return
      else
         open(file_unit,              &
              file   = file_name,     &
              status = 'new',         &
              form   = 'unformatted', &
              access = 'sequential')
         rewind(file_unit)
      end if

      write(file_unit) o_matrix

      close(file_unit, status = 'keep')

      if (allocated(o_matrix)) deallocate(o_matrix)

   end subroutine

   subroutine o_matrix_from_disc()

!     --------------------------------------------------------------------------
      integer :: p
      real(8) :: trace, diff, diff_tolerance
!     --------------------------------------------------------------------------

      inquire(file = file_name, exist = file_exists)

      if (file_exists) then
         open(file_unit,              &
              file   = file_name,     &
              status = 'old',         &
              form   = 'unformatted', &
              access = 'sequential')
         rewind(file_unit)
      else
!        panic
!        close(file_unit, status = 'keep')
!        return
      end if
      allocate(o_matrix(nr_mo_gerade + nr_mo_ungerade, nr_mo_gerade + nr_mo_ungerade))
      read(file_unit) o_matrix

      close(file_unit, status = 'keep')

      write(*, *) 'O_ia matrix (for PBHT MO Overlap Diagnostic)'
      call prqmat(o_matrix,       &
                  nr_mo,          &
                  nr_mo,          &
                  nr_mo,          &
                  nr_mo,          &
                  1,              &
                  (/1, 2, 3, 4/), &
                  6)
     
     trace = 0.0d0
     do p = 1, nr_mo
        trace = trace + o_matrix(p, p)
     end do
     write(*, *) 'trace(should be nr of orbitals):', trace
     write(*, *) 'trace/nr of orbitals (should be very close to 1):', trace/nr_mo
     diff = dabs((trace/nr_mo) - 1.0d0)
     diff_tolerance = 1.0d-4
     if (diff > diff_tolerance) then
        call quit('overlap_diagnostic: avg O_ia trace significantly different from 1.0')
     end if
   end subroutine

   subroutine integrate_o_matrix(block_length, w, ao)

!     --------------------------------------------------------------------------
      integer, intent(in)    :: block_length
      real(8), intent(in)    :: w(*)
      real(8), intent(in)    :: ao(block_length, *)
!     --------------------------------------------------------------------------
      integer                :: s, t
      integer                :: i, l, n, k, ic
      integer                :: iz, iq
      integer, external      :: iqfrompq
      real(8)                :: q_norm
      real(8),dimension(:,:),allocatable::q

!     --------------------------------------------------------------------------
!
!     structure of mo_coef (assuming no active orbitals):
!
!      _____  1 ... nr_quaternion_blocks
!     |  ___
!     | |  _  1 ... nr gerade mo negative energy secondary
!     | | |_  1 ... nr gerade ao (order according to ao_off)
!     | |___
!     |  ___
!     | |  _  1 ... nr gerade mo positive energy inactive
!     | | |_  1 ... nr gerade ao (order according to ao_off)
!     | |___
!     |  ___
!     | |  _  1 ... nr gerade mo positive energy secondary
!     | | |_  1 ... nr gerade ao (order according to ao_off)
!     | |___
!     |_____
!
!      _____  1 ... nr_quaternion_blocks
!     |  ___
!     | |  _  1 ... nr ungerade mo negative energy secondary
!     | | |_  1 ... nr ungerade ao (order according to ao_off)
!     | |___
!     |  ___
!     | |  _  1 ... nr ungerade mo positive energy inactive
!     | | |_  1 ... nr ungerade ao (order according to ao_off)
!     | |___
!     |  ___
!     | |  _  1 ... nr ungerade mo positive energy secondary
!     | | |_  1 ... nr ungerade ao (order according to ao_off)
!     | |___
!     |_____

      spinor = 0.0d0
!     construct the 4*nc parts of gerade spinor s in point k
      i = 0
      do iz = 1, nr_quaternion_blocks
!        gerade
         do s = 1, nr_mo_gerade
            do ic = 1,nc
              do n = 1, nr_g_ao_blocks(ic)
                iq = iqfrompq(n,1,ic,iz)
                do l = g_ao_start(n,ic), g_ao_end(n,ic)
                   i = i + 1
                   do k = 1, block_length
                      spinor(k, s, iq,ic) = spinor(k, s, iq, ic) + ao(k, l)*mo_coef(i)
                   end do
                end do
              end do
            enddo
         enddo
      enddo

!     construct the 4*nc parts of ungerade spinor s in point k
      if (nr_fermion_coreps > 1) then
         do iz = 1, nr_quaternion_blocks
!           ungerade
           do s = nr_mo_gerade + 1, nr_mo
               do ic = 1,nc
                 do n = 1, nr_u_ao_blocks(ic)
                   iq = iqfrompq(n,2,ic,iz)
                   do l = u_ao_start(n,ic), u_ao_end(n,ic)
                      i = i + 1
                      do k = 1, block_length
                        spinor(k, s, iq, ic) = spinor(k, s, iq, ic) + ao(k, l)*mo_coef(i) 
                      end do
                   end do
                 end do
               end do
           end do
         end do
      end if


     allocate(q(block_length,nr_mo))
      q=0.0d0 
!     construct the modulus of spinor s in point k
      do ic = 1,nc
        do iq = 1,4 !nr_quaternion_blocks 
          do s=1,nr_mo
            do k=1,block_length
               q(k, s)= q(k, s)+ spinor(k, s, iq, ic)*spinor(k, s, iq, ic)
            enddo
          enddo
        enddo
      enddo
      do s=1,nr_mo
        do k=1,block_length
          q(k,s)=dsqrt(q(k,s))
        enddo        
      enddo
!     combine moduli of spinors s and t to overlap diagnostic matrix
      do s = 1, nr_mo
         do t = s, nr_mo
            q_norm = 0.0d0
            do k = 1, block_length
              q_norm = q_norm + w(k)*q(k, s)*q(k, t)
            enddo
            o_matrix(s, t) = o_matrix(s, t) + q_norm
         end do
      end do
   deallocate(q)

   end subroutine

end module
