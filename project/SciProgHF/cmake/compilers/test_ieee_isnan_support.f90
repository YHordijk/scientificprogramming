! test program to see whether it is possible to use the ieee exception handling modules
program test
   use, intrinsic :: ieee_arithmetic

   if (ieee_support_nan()) then
      print *, 'ieee_arithmetic supported by compiler'
   endif
end program
