! test program to see whether it is possible to use the ieee exception handling modules
program test
   real (kind=kind(0.d0)) :: x
   x = 1.
   if (.not.isnan(x)) then
      print *, 'isnan intrinsic supported by compiler'
   endif
end program
