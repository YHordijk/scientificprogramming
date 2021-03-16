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

module interface_functional_read

   use dft_cfg
   use xc_derv
   use xcfun_1_0

   implicit none

   public parse_functional
   public set_xc_fun_alda
   public set_xc_fun_xalda
   public set_kin_fun_alda
   public set_kin_fun_xalda
   public consistency_after_dft_input
   public report_after_dft_input
   public set_hf_exchange_factor
 
   real(8), public :: hf_exchange_factor
   logical, public :: pure_hf_run
   logical, save   :: already_checked = .true.
   private

contains

   subroutine set_hf_exchange_factor(h)
      real(8), intent(in) :: h
      hf_exchange_factor = h
   end subroutine

  function get_weight(word, functional)

!   ----------------------------------------------------------------------------
    character(*), intent(in)  :: word
    character(*), intent(in)  :: functional
    real(8)                   :: get_weight
!   ----------------------------------------------------------------------------
    integer                   :: i
!   ----------------------------------------------------------------------------

    get_weight = 0.0d0

    if (word_contains(word, '=')) then

      i = index(word, '=')

      if (word(1:i-1) == functional) then
        read(word(i+1:len(word)), *) get_weight
      end if

    else

      if (word == functional) then
        get_weight = 1.0d0
      end if

    end if

  end function

   subroutine parse_functional(line,     &
                               f,        &
                               hfx_out,  &
                               mu_out,   &
                               beta_out, &
                               set_hf_exchange_factor)
   
!     radovan: line is made lowercase, then cut into words
!              if word contains "=" the weight is read
!              otherwise it is set to 1.0
!     
!              equivalent functionals:
!              b3lyp
!              slaterx=0.8 bcorrx=0.72 hfx=0.2 vwnc=0.19 lypc=0.81
      
!     --------------------------------------------------------------------------
      character(80)              :: line
      type(functional)           :: f
      real(8), intent(out)       :: hfx_out
      real(8), intent(out)       :: mu_out
      real(8), intent(out)       :: beta_out
      logical, intent(in)        :: set_hf_exchange_factor
!     --------------------------------------------------------------------------
      character(80), allocatable :: word_array(:)
      integer                    :: i, nr_words
      real(8)                    :: w, h
      integer                    :: u
!     --------------------------------------------------------------------------
      
      write(*, '(a)') '* Using the automatic differentiation xc functional:'
      write(*, '(a)') '   (weight: functional)'

      hfx_out  = 0.0d0
      mu_out   = 0.0d0
      beta_out = 0.0d0

      f%w = 0.0d0
      
      line = lowercase(line)
      
      nr_words = word_count(line)
      
      allocate(word_array(nr_words))
      read(line, *) (word_array(i), i = 1, nr_words)
      
      do i = 1, nr_words
      
      
!       composite functionals
!       =====================
      
        w = get_weight(word_array(i), 'lda')
        if (dabs(w) > tiny(0.0d0)) then
           f%w(XC_SLATERX) = f%w(XC_SLATERX) + w
           f%w(XC_VWN5C  ) = f%w(XC_VWN5C  ) + w
        end if
      
        w = get_weight(word_array(i), 'blyp')
        if (dabs(w) > tiny(0.0d0)) then
           f%w(XC_BECKEX) = f%w(XC_BECKEX) + w
           f%w(XC_LYPC  ) = f%w(XC_LYPC  ) + w
        end if
      
        w = get_weight(word_array(i), 'b3lyp')
        if (dabs(w) > tiny(0.0d0)) then
           f%w(XC_SLATERX   ) = f%w(XC_SLATERX   ) + 0.80d0*w
           f%w(XC_BECKECORRX) = f%w(XC_BECKECORRX) + 0.72d0*w
           f%w(XC_VWN5C     ) = f%w(XC_VWN5C     ) + 0.19d0*w
           f%w(XC_LYPC      ) = f%w(XC_LYPC      ) + 0.81d0*w
           hfx_out = hfx_out + 0.2d0*w
        end if
      
        w = get_weight(word_array(i), 'bhandh')
        if (dabs(w) > tiny(0.0d0)) then
           f%w(XC_SLATERX   ) = f%w(XC_SLATERX   ) + 0.5d0*w
           f%w(XC_LYPC      ) = f%w(XC_LYPC      ) +       w
           hfx_out = hfx_out + 0.5d0*w
        end if

        w = get_weight(word_array(i), 'pp86')
        if (dabs(w) > tiny(0.0d0)) then
           f%w(XC_PW86X) = f%w(XC_PW86X) + w
           f%w(XC_P86C)  = f%w(XC_P86C)  + w
        end if

        w = get_weight(word_array(i), 'm05')
        if (dabs(w) > tiny(0.0d0)) then
           h = 0.28d0
!radovan:  notice the difference of the x weight in m05 and m06
!          it seems that the hfx modification is already absorbed
!          in the parameters of m06x
!          this is not the case for m05x
           f%w(XC_M05X) = f%w(XC_M05X) + w*(1.0d0 - h)
           f%w(XC_M05C) = f%w(XC_M05C) + w
           hfx_out = hfx_out + h*w
        end if
      
        w = get_weight(word_array(i), 'm06')
        if (dabs(w) > tiny(0.0d0)) then
           h = 0.27d0
           f%w(XC_M06X) = f%w(XC_M06X) + w
           f%w(XC_M06C) = f%w(XC_M06C) + w
           hfx_out = hfx_out + h*w
#ifndef MOD_UNRELEASED
           call quit('this functional is not available in this version')
#endif
        end if
      
        w = get_weight(word_array(i), 'm06l')
        if (dabs(w) > tiny(0.0d0)) then
           f%w(XC_M06LX) = f%w(XC_M06LX) + w
           f%w(XC_M06LC) = f%w(XC_M06LC) + w
#ifndef MOD_UNRELEASED
           call quit('this functional is not available in this version')
#endif
        end if
      
        w = get_weight(word_array(i), 'm06-l')
        if (dabs(w) > tiny(0.0d0)) then
           f%w(XC_M06LX) = f%w(XC_M06LX) + w
           f%w(XC_M06LC) = f%w(XC_M06LC) + w
#ifndef MOD_UNRELEASED
           call quit('this functional is not available in this version')
#endif
        end if
      
        w = get_weight(word_array(i), 'pbe')
        if (dabs(w) > tiny(0.0d0)) then
          write(*, '(f10.6, a)') w, ': pbe'
           f%w(XC_PBEX) = f%w(XC_PBEX) + w
           f%w(XC_PBEC) = f%w(XC_PBEC) + w
        end if
      
        w = get_weight(word_array(i), 'pbe0')
        if (dabs(w) > tiny(0.0d0)) then
           f%w(XC_PBEX) = f%w(XC_PBEX) + 0.75d0*w
           f%w(XC_PBEC) = f%w(XC_PBEC) +        w
           hfx_out = hfx_out + 0.25d0*w
        end if
      
        w = get_weight(word_array(i), 'pbe38')
        if (dabs(w) > tiny(0.0d0)) then
           f%w(XC_PBEX) = f%w(XC_PBEX) + 0.625d0*w
           f%w(XC_PBEC) = f%w(XC_PBEC) +        w
           hfx_out = hfx_out + 0.375d0*w
        end if
      
        w = get_weight(word_array(i), 'kt1')
        if (dabs(w) > tiny(0.0d0)) then
           f%w(XC_SLATERX) = f%w(XC_SLATERX) +            w
           f%w(XC_KTX    ) = f%w(XC_KTX    ) -    0.006d0*w
           f%w(XC_VWN5C  ) = f%w(XC_VWN5C  ) +            w
        end if
      
        w = get_weight(word_array(i), 'kt2')
        if (dabs(w) > tiny(0.0d0)) then
           f%w(XC_SLATERX) = f%w(XC_SLATERX) +  1.07173d0*w
           f%w(XC_KTX    ) = f%w(XC_KTX    ) -    0.006d0*w
           f%w(XC_VWN5C  ) = f%w(XC_VWN5C  ) + 0.576727d0*w
        end if
      
      
!       "atomic" functionals
!       ====================

        w = get_weight(word_array(i), 'pw86x')
        if (dabs(w) > tiny(0.0d0)) then
           f%w(XC_PW86X) = f%w(XC_PW86X) + w
        end if

        w = get_weight(word_array(i), 'p86c')
        if (dabs(w) > tiny(0.0d0)) then
           f%w(XC_P86C) = f%w(XC_P86C) + w
        end if
      
        w = get_weight(word_array(i), 'kin_pw91')
        if (dabs(w) > tiny(0.0d0)) then
          write(*, '(f10.6, a)') w, ': kin_pw91'
          f%w(XC_PW91K) = f%w(XC_PW91K) + w
        end if

        w = get_weight(word_array(i), 'kin_tf')
        if (dabs(w) > tiny(0.0d0)) then
          write(*, '(f10.6, a)') w, ': kin_tf'
          f%w(XC_TFK) = f%w(XC_TFK) + w
        end if

        w = get_weight(word_array(i), 'ktx')
        if (dabs(w) > tiny(0.0d0)) then
           f%w(XC_KTX) = f%w(XC_KTX) + w
        end if
      
        w = get_weight(word_array(i), 'pbex')
        if (dabs(w) > tiny(0.0d0)) then
          write(*, '(f10.6, a)') w, ': pbex'
           f%w(XC_PBEX) = f%w(XC_PBEX) + w
        end if
      
        w = get_weight(word_array(i), 'pbec')
        if (dabs(w) > tiny(0.0d0)) then
          write(*, '(f10.6, a)') w, ': pbec'
           f%w(XC_PBEC) = f%w(XC_PBEC) + w
        end if
      
        w = get_weight(word_array(i), 'lypc')
        if (dabs(w) > tiny(0.0d0)) then
           f%w(XC_LYPC) = f%w(XC_LYPC) + w
        end if
      
        w = get_weight(word_array(i), 'hfx')
        if (dabs(w) > tiny(0.0d0)) then
           hfx_out = hfx_out + w
        end if
      
        w = get_weight(word_array(i), 'ldaerfx')
        if (dabs(w) > tiny(0.0d0)) then
           f%w(XC_LDAERFX) = f%w(XC_LDAERFX) + w
        end if
      
        w = get_weight(word_array(i), 'ldaerfc')
        if (dabs(w) > tiny(0.0d0)) then
           f%w(XC_LDAERFC) = f%w(XC_LDAERFC) + w
        end if
      
        w = get_weight(word_array(i), 'ldaerfc_jt')
        if (dabs(w) > tiny(0.0d0)) then
           f%w(XC_LDAERFC_JT) = f%w(XC_LDAERFC_JT) + w
        end if
      
        w = get_weight(word_array(i), 'mu')
        if (dabs(w) > tiny(0.0d0)) then
           f%w(XC_RANGESEP_MU) = w
           mu_out = mu_out + w
        end if
      
        w = get_weight(word_array(i), 'beta')
        if (dabs(w) > tiny(0.0d0)) then
           beta_out = beta_out + w
        end if

        w = get_weight(word_array(i), 'slaterx')
        if (dabs(w) > tiny(0.0d0)) then
           f%w(XC_SLATERX) = f%w(XC_SLATERX) + w
        end if

        w = get_weight(word_array(i), 'vwnc')
        if (dabs(w) > tiny(0.0d0)) then
           f%w(XC_VWN5C) = f%w(XC_VWN5C) + w
        end if

        w = get_weight(word_array(i), 'beckex')
        if (dabs(w) > tiny(0.0d0)) then
           f%w(XC_BECKEX) = f%w(XC_BECKEX) + w
        end if

        w = get_weight(word_array(i), 'bcorrx')
        if (dabs(w) > tiny(0.0d0)) then
           f%w(XC_BECKECORRX) = f%w(XC_BECKECORRX) + w
        end if
      
      end do
      
      deallocate(word_array)
 
      write(*, '(a)')

      if (set_hf_exchange_factor) then
         hf_exchange_factor = hfx_out
      end if

!     check that not all weights are zero
      h = 0.0d0
      do i = 1, size(f%w)
         h = h + dabs(f%w(i)) 
      end do
      if (.not.already_checked) pure_hf_run = .false.
      if (h < tiny(0.0d0)) then
         if (hf_exchange_factor == 1.0d0) then
             pure_hf_run = .true.
             already_checked = .true.
         else
!            no functional is recognized
!            in this case quit the program
             print *, 'XCFun functional not recognized from input'
             stop
         end if
      end if
   
   end subroutine

   subroutine set_xc_fun_alda(f)

!     --------------------------------------------------------------------------
      type(functional) :: f
!     --------------------------------------------------------------------------

      f%w = 0.0d0
      f%w(XC_SLATERX) = 1.0d0
      f%w(XC_VWN5C  ) = 1.0d0

   end subroutine

   subroutine set_xc_fun_xalda(f)

!     --------------------------------------------------------------------------
      type(functional) :: f
!     --------------------------------------------------------------------------

      f%w = 0.0d0
      f%w(XC_SLATERX) = 1.0d0 - hf_exchange_factor
      f%w(XC_VWN5C  ) = 1.0d0

   end subroutine

   subroutine set_kin_fun_alda(f)

!     --------------------------------------------------------------------------
      type(functional) :: f
!     --------------------------------------------------------------------------

      f%w = 0.0d0
      f%w(XC_TFK) = 1.0d0

   end subroutine

   subroutine set_kin_fun_xalda(f)

!     --------------------------------------------------------------------------
      type(functional) :: f
!     --------------------------------------------------------------------------

      call set_kin_fun_alda(f)
   end subroutine

  subroutine consistency_after_dft_input()

    if (dft_cfg_asymptote_is_lb94 .and. dft_cfg_asymptote_is_lbalpha) then
       print *, 'LBalpha = LB94 = .true.'
       stop
    end if

    if (dft_cfg_sdft_collinear .and. dft_cfg_no_sdft) then
       print *, '.COLLINEAR cannot be used together with .NOSDFT'
       stop
    end if

  end subroutine

  subroutine report_after_dft_input()

!   ----------------------------------------------------------------------------
    integer :: u
!   ----------------------------------------------------------------------------

    write(*, *) ' ===== Kohn-Sham calculation set-up ====='

!   thresholds

    write(*, '(a)')          ' * DFT thresholds:'
    write(*, '(a, e12.5)')   '   - small density threshold = ', dft_cfg_tinydens

!   asymptotic correction

    if (dft_cfg_grac) then
      write(*, '(a)')        ' * Gradient regulated asymptotic correction:'
      write(*, '(a, e12.5)') '   - alpha = ', dft_cfg_grac_alpha
      write(*, '(a, e12.5)') '   - beta  = ', dft_cfg_grac_beta
      write(*, '(a)')        '   - threshold for difference in HOMO eigenvalue'
      write(*, '(a, e12.5)') '     below which asymptotic correction '// &
                             'is switched on: ', dft_cfg_ac_threshold
      write(*, '(a, e12.5)') '   - bulk potential will be shifted by (IP): ', dft_cfg_ac_ip
    end if

    if (dft_cfg_saop) then
      write(*, '(a)')        ' * Statistical averaging of (model) orbital potentials:'
    end if

    if (dft_cfg_saop_with_response_part) then
      write(*, '(a)')        '   - original SAOP as in JCP 112, 1344 (2000) chosen'
      write(*, '(a)')        '   - this activates ALDA xc kernel'
      write(*, '(a)')        '   - the functional under .DFT is expected to be GLLBhole'
    end if

    if (dft_cfg_asymptote_is_lb94) then
      write(*, '(a)')        '   - the asymptotic potential is LB94'
    end if

    if (dft_cfg_asymptote_is_lbalpha) then
      write(*, '(a)')        '   - the asymptotic potential is LBalpha'
    end if

!   sdft

    if (dft_cfg_no_sdft) then

      write(*, '(a)')        ' * No spin density contribution in response calculations.'

    else

      write(*, '(a)')        ' * Spin density contribution in response calculations:'

      if (dft_cfg_sdft_collinear) then
        write(*, '(a)')      '   - use the collinear approximation as '// &
                             'a definition of the spin density'
      else
        write(*, '(a)')      '   - use the norm of the spin magnetization vector as '// &
                             'a definition of the spin density (noncollinear definition)'
      end if

    end if !if (dft_cfg_no_sdft) then

!   alda

    if (dft_cfg_alda_hs .or. dft_cfg_alda_ha) then

      write(*, '(a)')        ' * Adiabatic local density approximation (ALDA):'
      write(*, '(a)')        '   - approximate all functional derivatives '// &
                             'beyond the xc potential by SVWN derivatives'

      if (dft_cfg_alda_hs .and. .not. dft_cfg_alda_ha) then
        write(*, '(a)')      '   - use ALDA only for the Hermitian part '// &
                             '(density contribution)'
        write(*, '(a)')      '   - use the proper xc kernel for the anti-Hermitian '// &
                             'part (spin density contribution)'
      end if

      if (dft_cfg_alda_ha .and. .not. dft_cfg_alda_hs) then
        write(*, '(a)')      '   - use ALDA only for the anti-Hermitian part '// &
                             '(spin density contribution)'
        write(*, '(a)')      '   - use the proper xc kernel for the Hermitian '// &
                             'part (density contribution)'
      end if

      if (dft_cfg_xalda) then
        write(*, '(a)')      '   - keep the fraction of exact exchange of the '// &
                             'xc functional in the solution of the response equation'
      else
        write(*, '(a)')      '   - for hybrid functionals exact exchange is '// &
                             'switched off in the solution of the response equation'
        write(*, '(a)')      '     (if you want to keep it, use .XALDA)'
      end if

    end if !if (dft_cfg_alda_hs .or. dft_cfg_alda_ha) then

  end subroutine

  function word_count(s)

!   ----------------------------------------------------------------------------
    character(*), intent(in) :: s
    integer                  :: word_count
!   ----------------------------------------------------------------------------
    integer                  :: i
    logical                  :: is_blank
!   ----------------------------------------------------------------------------

    word_count = 0

    if (len(s) <= 0) return

    is_blank = .true.

    do i = 1, len(s)
      if (s(i:i) == ' ') then
        is_blank = .true.
      else if (is_blank) then
        word_count = word_count + 1
        is_blank = .false.
      end if
    end do

  end function

  function word_contains(word, substring)

!   ----------------------------------------------------------------------------
    character(*), intent(in) :: word
    character(*), intent(in) :: substring
!   ----------------------------------------------------------------------------
    logical                  :: word_contains
!   ----------------------------------------------------------------------------

    word_contains = .false.
    if (index(word, substring) > 0) then
      word_contains = .true.
    end if

  end function

  function lowercase(s)

!   ----------------------------------------------------------------------------
    character(*), intent(in) :: s
    character(len(s))        :: lowercase
!   ----------------------------------------------------------------------------
    integer                  :: off, i, ia
!   ----------------------------------------------------------------------------

    lowercase = s

    off = iachar('a') - iachar('A')

    do i = 1, len(s)
      ia = iachar(s(i:i))
      if (ia >= iachar('A') .and. ia <= iachar('Z')) then
        lowercase(i:i) = achar(ia + off)
      end if
    enddo

  end function

end module
