program i8_test
!
! This program ALWAYS crashes when linked against integer*8 blas, but passes with integer*4 blas
!
! Invented by Hans-Jorgen Aa. Jensen, SDU Odense, Denmark
!
! Adapted into Dirac CMake system by Miro Ilias, GSI, Darmstadt
!
  integer*8 N_8
  integer*4 N_4(2)
  equivalence (N_8, N_4)
  real*8    a(20)
  integer*8 j8, idamax
  a(:) = 0.0d0
  a(5) = 1.0d0
  N_4(1) = 20
  N_4(2) = 20
  j8 = idamax(N_8,a,1)
end
