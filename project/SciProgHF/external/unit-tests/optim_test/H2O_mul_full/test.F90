program test

  use second_order_minimization
  use unit_testing
  implicit none

  call minimization_drv( unit_test = .true. )

  call compare_files(60,'OPTIM.OUT','OPTIM.REF',1.0d-9)

  if(all_tests_passed()) write(6,*) ' test result: =====  OK  ====='

end program
