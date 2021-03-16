module machine_parameters

      implicit none
      real(8), external :: dlamch
      real(8), public :: mp_eps   
      real(8), public :: mp_sfmin 
      real(8), public :: mp_base  
      real(8), public :: mp_prec  
      real(8), public :: mp_t     
      real(8), public :: mp_rnd   
      real(8), public :: mp_emin  
      real(8), public :: mp_rmin  
      real(8), public :: mp_emax  
      real(8), public :: mp_rmax  

contains

      subroutine get_machine_parameters

      real(8), external :: dlamch

      mp_eps   = DLAMCH( 'Epsilon' )
      mp_sfmin = DLAMCH( 'Safe minimum' )
      mp_base  = DLAMCH( 'Base' )
      mp_prec  = DLAMCH( 'Precision' )
      mp_t     = DLAMCH( 'Number of digits in mantissa' )
      mp_rnd   = DLAMCH( 'Rounding mode' )
      mp_emin  = DLAMCH( 'Minimum exponent' )
      mp_rmin  = DLAMCH( 'Underflow threshold' )
      mp_emax  = DLAMCH( 'Largest exponent' )
      mp_rmax  = DLAMCH( 'Overflow threshold' )

      end subroutine get_machine_parameters

      subroutine print_machine_parameters

      WRITE( 6, * )' Epsilon                      = ', mp_eps
      WRITE( 6, * )' Safe minimum                 = ', mp_sfmin
      WRITE( 6, * )' Base                         = ', mp_base
      WRITE( 6, * )' Precision                    = ', mp_prec
      WRITE( 6, * )' Number of digits in mantissa = ', mp_t
      WRITE( 6, * )' Rounding mode                = ', mp_rnd
      WRITE( 6, * )' Minimum exponent             = ', mp_emin
      WRITE( 6, * )' Underflow threshold          = ', mp_rmin
      WRITE( 6, * )' Largest exponent             = ', mp_emax
      WRITE( 6, * )' Overflow threshold           = ', mp_rmax
      WRITE( 6, * )' Reciprocal of safe minimum   = ', 1 / mp_sfmin

      end subroutine print_machine_parameters

end module
