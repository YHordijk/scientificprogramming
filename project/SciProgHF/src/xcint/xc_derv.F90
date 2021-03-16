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

module xc_derv

   use xc_max_block_length
   use xc_mpi
   use xcfun_1_0
   use dft_cfg
 
   implicit none
 
   public get_functional_derv
   public fun_is_lda
   public fun_is_gga
   public fun_is_tau_mgga
   public xcfun_set_functional

   public functional
   type functional
      real(8)              :: w(XC_NR_PARAMS)
      integer              :: id
      logical              :: xcfun_id_is_set
      real(8), allocatable :: input(:, :)
      real(8), allocatable :: output(:, :)
   end type

   logical, public :: fun_is_automatic

   type(functional), public :: xc_fun
   type(functional), public :: xc_fun_alda
   type(functional), public :: xc_fun_xalda

   integer, parameter, public :: d0000000 =  0
   integer, parameter, public :: d1000000 =  1
   integer, parameter, public :: d0010000 =  2
   integer, parameter, public :: d0000100 =  3
   integer, parameter, public :: d2000000 =  4
   integer, parameter, public :: d0200000 =  5
   integer, parameter, public :: d0020000 =  6
   integer, parameter, public :: d0002000 =  7
   integer, parameter, public :: d1010000 =  8
   integer, parameter, public :: d0101000 =  9
   integer, parameter, public :: d1000100 = 10
   integer, parameter, public :: d0010100 = 11
   integer, parameter, public :: d3000000 = 12
   integer, parameter, public :: d0030000 = 13
   integer, parameter, public :: d1200000 = 14
   integer, parameter, public :: d1020000 = 15
   integer, parameter, public :: d1002000 = 16
   integer, parameter, public :: d2010000 = 17
   integer, parameter, public :: d0210000 = 18
   integer, parameter, public :: d1101000 = 19
   integer, parameter, public :: d0111000 = 20
   integer, parameter, public :: d0012000 = 21
   
   integer, parameter, public :: d0000010 = 22
   integer, parameter, public :: d0000020 = 23
   integer, parameter, public :: d1000010 = 24
   integer, parameter, public :: d0010010 = 25
   integer, parameter, public :: d0100001 = 26
   integer, parameter, public :: d0001001 = 27
   integer, parameter, public :: d0000002 = 28
   
   integer, parameter, public :: nr_nonzero_derv = 28
   
   private

#ifdef PRG_DIRAC
#include "../cfun/functionals.h"
#endif

   integer, parameter :: max_dim = 99

contains

   function fun_is_lda(f)
      type(functional) :: f
      logical :: fun_is_lda

      fun_is_lda = .false.
      if (fun_is_automatic) then
         fun_is_lda = (xc_get_type(f%id) == 0)
      else
#ifdef PRG_DIRAC
         fun_is_lda = (.not. dftisgga())
#endif
      end if
   end function

   function fun_is_gga(f)
      type(functional) :: f
      logical :: fun_is_gga

      fun_is_gga = .false.
      if (fun_is_automatic) then
         fun_is_gga = (xc_get_type(f%id) == 1)
      else
#ifdef PRG_DIRAC
         fun_is_gga = dftisgga()
#endif
      end if
   end function

   function fun_is_tau_mgga(f)
      type(functional) :: f
      logical :: fun_is_tau_mgga

      fun_is_tau_mgga = .false.
      if (fun_is_automatic) then
         fun_is_tau_mgga = (xc_get_type(f%id) == 2)
      else
         fun_is_tau_mgga = .false.
      end if
   end function

   subroutine xcfun_set_functional(f,              &
                                   response_order, &
                                   parallel_xc)

!     --------------------------------------------------------------------------
      type(functional), intent(inout) :: f
      integer,          intent(in)    :: response_order
      logical,          intent(in)    :: parallel_xc
!     --------------------------------------------------------------------------
      integer                         :: i, l
!     --------------------------------------------------------------------------

      if (.not. f%xcfun_id_is_set) then
#ifdef VAR_MPI
         if (parallel_xc) then
            call xc_mpi_bcast(f%w)
            call xc_mpi_bcast(f%id)
            call xc_mpi_bcast(f%xcfun_id_is_set)
         end if
#endif
         f%id = xc_new_functional()
         call xc_set_mode(f%id, XC_VARS_NS)

         do i = 1, XC_NR_PARAMS
            if (dabs(f%w(i)) > tiny(0.0d0)) then
               call xc_set_param(f%id, i, f%w(i))
            end if
         end do
         f%xcfun_id_is_set = .true.
      end if

      if (allocated(f%input)) then
         l = size(f%input)/max_block_length
      else
         l = 0
      end if
      if (xc_input_length(f%id) > l) then
         l = xc_input_length(f%id)
         if (allocated(f%input)) then
            deallocate(f%input)
         end if
         allocate(f%input(l, max_block_length))
      end if

      if (allocated(f%output)) then
         l = size(f%output)/max_block_length
      else
         l = 0
      end if
      if (xc_output_length(f%id, 1 + response_order) > l) then
         l = xc_output_length(f%id, 1 + response_order)
         if (allocated(f%output)) then
            deallocate(f%output)
         end if
         allocate(f%output(l, max_block_length))
      end if

   end subroutine

  subroutine get_functional_derv(f,            &
                                 f_alda,       &
                                 f_xalda,      &
                                 order,        &
                                 block_length, &
                                 w,            &
                                 r_0,          &
                                 z_0,          &
                                 n,            &
                                 derv,         &
                                 alda_real,    &
                                 alda_imag,    &
                                 xalda,        &
                                 tau)

!   ----------------------------------------------------------------------------
    type(functional)              :: f
    type(functional)              :: f_alda, f_xalda
    integer,          intent(in)  :: order
    integer,          intent(in)  :: block_length
    real(8),          intent(in)  :: w(*)
    real(8),          intent(in)  :: r_0(*)
    real(8),          intent(in)  :: z_0(*)
    integer,          intent(in)  :: n
    real(8),          intent(out) :: derv(max_block_length, 0:n)
    logical, optional, intent(in) :: alda_real
    logical, optional, intent(in) :: alda_imag
    logical, optional, intent(in) :: xalda
    real(8), optional, intent(in) :: tau(*)
!   ----------------------------------------------------------------------------
    real(8)                       :: t1(max_dim, max_block_length)
    real(8)                       :: t2(max_dim, max_block_length)
    integer                       :: k, l, id
#ifdef PRG_DIRAC
    real(8), external             :: dftenergy
#endif
    logical                       :: do_alda_real
    logical                       :: do_alda_imag
    logical                       :: do_xalda
!   ----------------------------------------------------------------------------

    do_alda_real = .false.
    do_alda_imag = .false.
    do_xalda     = .false.

    if (present(alda_real)) do_alda_real = alda_real
    if (present(alda_imag)) do_alda_imag = alda_imag
    if (present(xalda))     do_xalda     = xalda

!   if less then linear response ignore alda
    if (order < 2) then
       do_alda_real = .false.
       do_alda_imag = .false.
    end if

    derv = 0.0d0
    t1   = 0.0d0
    t2   = 0.0d0

    if (fun_is_automatic) then

      f%input = 0.0d0
      do k = 1, block_length
         f%input(1, k) = r_0(k)
      end do
      if (fun_is_gga(f)) then
         do k = 1, block_length
            f%input(3, k) = z_0(k)
         end do
      end if
      if (fun_is_tau_mgga(f)) then
         do k = 1, block_length
            f%input(3, k) = z_0(k)
         end do
         do k = 1, block_length
            f%input(6, k) = tau(k)
         end do
      end if

      call xc_eval(f%id,         &
                   order,        &
                   block_length, &
                   f%input,      &
                   f%output)

      if (do_alda_real .or. do_alda_imag) then
         f_alda%input = 0.0d0
         do k = 1, block_length
            f_alda%input(1, k) = r_0(k)
         end do
         id = f_alda%id
         if (do_xalda) then
            id = f_xalda%id
         end if
         call xc_eval(id,           &
                      order,        &
                      block_length, &
                      f_alda%input, &
                      f_alda%output)
      end if

      do k = 1, block_length

         derv(k, d0000000) = f%output(1, k)

         if (order == 0) then
               call my_quit('get_functional_derv: order 0; probably programming error')
         end if

         if (order > 0) then
               if (fun_is_lda(f)) then
                  derv(k, d1000000) = f%output(XC_D10, k)
               else if (fun_is_gga(f)) then
                  derv(k, d1000000) = f%output(XC_D10000, k)
                  derv(k, d0010000) = f%output(XC_D00100, k)
               else if (fun_is_tau_mgga(f)) then
                  derv(k, d1000000) = f%output(XC_D1000000, k)
                  derv(k, d0010000) = f%output(XC_D0010000, k)
                  derv(k, d0000010) = f%output(XC_D0000010, k)
               end if
         end if

         if (order > 1) then
!              order 2 - spin unpolarized
               if (do_alda_real) then
                  derv(k, d2000000) = f_alda%output(XC_D20, k)
               else
                  if (fun_is_lda(f)) then
                     derv(k, d2000000) = f%output(XC_D20, k)
                  else if (fun_is_gga(f)) then
                     derv(k, d2000000) = f%output(XC_D20000, k)
                     derv(k, d1010000) = f%output(XC_D10100, k)
                     derv(k, d0020000) = f%output(XC_D00200, k)
                     derv(k, d0010000) = f%output(XC_D00100, k)
                  else if (fun_is_tau_mgga(f)) then
                     derv(k, d0010000) = f%output(XC_D0010000, k)
                     derv(k, d2000000) = f%output(XC_D2000000, k)
                     derv(k, d1010000) = f%output(XC_D1010000, k)
                     derv(k, d0020000) = f%output(XC_D0020000, k)

                     derv(k, d0000020) = f%output(XC_D0000020, k)
                     derv(k, d0010010) = f%output(XC_D0010010, k)
                     derv(k, d1000010) = f%output(XC_D1000010, k)
                  end if
               end if
!              order 2 - spin polarized
               if (do_alda_imag) then
                  derv(k, d0200000) = f_alda%output(XC_D02, k)
               else
                  if (fun_is_lda(f)) then
                     derv(k, d0200000) = f%output(XC_D02, k)
                  else if (fun_is_gga(f)) then
                     derv(k, d0200000) = f%output(XC_D02000, k)
                     derv(k, d0101000) = f%output(XC_D01010, k)
                     derv(k, d0002000) = f%output(XC_D00020, k)
                     derv(k, d0000100) = f%output(XC_D00001, k)
                  else if (fun_is_tau_mgga(f)) then
                     derv(k, d0000100) = f%output(XC_D0000100, k)
                     derv(k, d0200000) = f%output(XC_D0200000, k)
                     derv(k, d0101000) = f%output(XC_D0101000, k)
                     derv(k, d0002000) = f%output(XC_D0002000, k)

                     derv(k, d0000002) = f%output(XC_D0000002, k)
                     derv(k, d0001001) = f%output(XC_D0001001, k)
                     derv(k, d0100001) = f%output(XC_D0100001, k)
                  end if
               end if
         end if

         if (order > 2) then
               if (do_alda_real .or. do_alda_imag) then
                  call my_quit('get_functional_derv: xcfun alda qr not implemented')
               end if
               if (fun_is_lda(f)) then
                  derv(k, d2000000) = f%output(XC_D20, k)
                  derv(k, d0200000) = f%output(XC_D02, k)
                  derv(k, d1200000) = f%output(XC_D12, k)
                  derv(k, d3000000) = f%output(XC_D30, k)
               else if (fun_is_gga(f)) then
                  derv(k, d0010000) = f%output(XC_D00100, k)
                  derv(k, d0000100) = f%output(XC_D00001, k)
                  derv(k, d2000000) = f%output(XC_D20000, k)
                  derv(k, d1010000) = f%output(XC_D10100, k)
                  derv(k, d1000100) = f%output(XC_D10001, k)
                  derv(k, d0200000) = f%output(XC_D02000, k)
                  derv(k, d0101000) = f%output(XC_D01010, k)
                  derv(k, d0020000) = f%output(XC_D00200, k)
                  derv(k, d0010100) = f%output(XC_D00101, k)
                  derv(k, d0002000) = f%output(XC_D00020, k)
                  derv(k, d3000000) = f%output(XC_D30000, k)
                  derv(k, d2010000) = f%output(XC_D20100, k)
                  derv(k, d1200000) = f%output(XC_D12000, k)
                  derv(k, d1101000) = f%output(XC_D11010, k)
                  derv(k, d1020000) = f%output(XC_D10200, k)
                  derv(k, d1002000) = f%output(XC_D10020, k)
                  derv(k, d0210000) = f%output(XC_D02100, k)
                  derv(k, d0111000) = f%output(XC_D01110, k)
                  derv(k, d0030000) = f%output(XC_D00300, k)
                  derv(k, d0012000) = f%output(XC_D00120, k)
               else if (fun_is_tau_mgga(f)) then
                  call my_quit('get_functional_derv: tau mgga qr not implemented')
               end if
         end if

         if (order > 3) then
               call my_quit('get_functional_derv: order too high')
         end if

      end do

    else

#ifdef PRG_DIRAC
       do k = 1, block_length
          if (r_0(k) > dft_cfg_tinydens) then
             select case (order)
                case (0)
                   derv(k, d0000000) = dftenergy(r_0(k), dsqrt(z_0(k)))
                case (1)
                   derv(k, d0000000) = dftenergy(r_0(k), dsqrt(z_0(k)))
                   call dftpot0(t1(1, k), 1.0d0, r_0(k), dsqrt(z_0(k)))
                case (2)
                   call  dftpot1(t1(1, k), 1.0d0, r_0(k), dsqrt(z_0(k)))
                   call sdftpot1(t2(1, k), 1.0d0, r_0(k), dsqrt(z_0(k)))
                case (3)
                   call sdftpot2(t1(1, k), 1.0d0, r_0(k), dsqrt(z_0(k)))
                case default
                   call my_quit('get_functional_derv: order too high')
             end select
          end if
       end do
#else
      print *, 'tried to use cfun in standalone'
      stop
#endif

      do l = 1, n
        do k = 1, block_length
          derv(k, l) = derv(k, l) + t1(l, k)
        end do
      end do
      if (order == 2) then
         do l = 1, n
           do k = 1, block_length
             derv(k, l) = derv(k, l) + t2(l, k)
           end do
         end do
      end if

    end if

!   multiply by weight
    do l = 0, n
       do k = 1, block_length
          derv(k, l) = w(k)*derv(k, l)
       end do
    end do

  end subroutine


end module
