

subroutine check_for_nan(var,caller)
#ifdef HAVE_IEEE_ISNAN
   use, intrinsic :: ieee_arithmetic
#endif
   character(6)   :: caller 
   real*8         :: var

#ifdef HAVE_IEEE_ISNAN
   if (ieee_is_nan(var)) then
#elif HAVE_INTRINSIC_ISNAN
   if (isnan(var)) then
#else
   if (.false.) then
#endif
      call quit("Stopping, NaN detected at "//caller(1:6))
   endif

end subroutine

