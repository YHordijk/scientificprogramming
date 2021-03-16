module qstack
#if  defined (INT_STAR8)
#define  CINT48 c_long_long
#else
#define  CINT48 c_int
#endif

       interface

         FUNCTION qstack_init (i1,i2) bind(c, name="qstack_init")
           use, intrinsic ::  iso_c_binding
           integer (CINT48)         :: qstack_init
           integer (CINT48)         :: i1,i2
         END FUNCTION
       
         FUNCTION qstack_push (i1, i2, ra1) bind(c, name="qstack_push")
           use, intrinsic ::  iso_c_binding
           integer (CINT48)         :: qstack_push
           integer (CINT48)         :: i1, i2
           real (c_double)          :: ra1(*)
         END FUNCTION
       
         FUNCTION qstack_popf (i1, ra1) bind(c, name="qstack_popf")
           use, intrinsic ::  iso_c_binding
           integer (CINT48)         :: qstack_popf
           integer (CINT48)         :: i1
           real (c_double)          :: ra1(*)
         END FUNCTION
       
         FUNCTION qstack_popl (i1, ra1) bind(c, name="qstack_popl")
           use, intrinsic ::  iso_c_binding
           integer (CINT48)         :: qstack_popl
           integer (CINT48)         :: i1
           real (c_double)          :: ra1(*)
         END FUNCTION
       
         FUNCTION qstack_peekf (i1, ra1) bind(c, name="qstack_peekf")
           use, intrinsic ::  iso_c_binding
           integer (CINT48)         :: qstack_peekf
           integer (CINT48)         :: i1
           real (c_double)          :: ra1(*)
         END FUNCTION
       
         FUNCTION qstack_peekl (i1, ra1) bind(c, name="qstack_peekl")
           use, intrinsic ::  iso_c_binding
           integer (CINT48)         :: qstack_peekl
           integer (CINT48)         :: i1
           real (c_double)          :: ra1(*)
         END FUNCTION
       
         FUNCTION qstack_peekn (i1, ra1, i2) bind(c, name="qstack_peekn")
           use, intrinsic ::  iso_c_binding
           integer (CINT48)         :: qstack_peekn
           integer (CINT48)         :: i1
           real (c_double)          :: ra1(*)
           integer (CINT48)         :: i2
         END FUNCTION
       
         FUNCTION qstack_compl (i1, ra1) bind(c, name="qstack_compl")
           use, intrinsic ::  iso_c_binding
           integer (CINT48)         :: qstack_compl
           integer (CINT48)         :: i1
           real (c_double)          :: ra1(*)
         END FUNCTION
       
         FUNCTION qstack_getinf () bind(c, name="qstack_getinf")
           use, intrinsic ::  iso_c_binding
           integer (CINT48)         :: qstack_getinf
         END FUNCTION
       
         FUNCTION qstack_meminf () bind(c, name="qstack_meminf")
           use, intrinsic ::  iso_c_binding
           real (c_double)         :: qstack_meminf
         END FUNCTION
       
         FUNCTION qstack_drop (i1) bind(c, name="qstack_drop")
           use, intrinsic ::  iso_c_binding
           integer (CINT48)         :: qstack_drop
           integer (CINT48)         :: i1
         END FUNCTION
       
         FUNCTION qstack_directx_r (i1, i2, ra1, i3, r1) bind(c, name="qstack_directx_r")
           use, intrinsic ::  iso_c_binding
           integer (CINT48)         :: qstack_directx_r
           integer (CINT48)         :: i1
           integer (CINT48)         :: i2
           real (c_double)          :: ra1(*)
           integer (CINT48)         :: i3
           real (c_double)          :: r1
         END FUNCTION
       
         FUNCTION qstack_directx_c (i1, i2, ra1, i3, c1) bind(c, name="qstack_directx_c")
           use, intrinsic ::  iso_c_binding
           integer (CINT48)                :: qstack_directx_c
           integer (CINT48)                :: i1
           integer (CINT48)                :: i2
           real (c_double)                 :: ra1(*)
           integer (CINT48)                :: i3
           complex (c_double_complex)      :: c1
         END FUNCTION
       
         SUBROUTINE qstack_shutdown () bind(c, name="qstack_shutdown")
           use, intrinsic :: iso_c_binding
         END SUBROUTINE

       end interface

end module
