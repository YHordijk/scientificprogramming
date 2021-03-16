subroutine xcint_potential_rks(mat_dim, &
                               dmat_0,  &
                               fmat)

  use xcint_main
 ! use include_dgroup_h,only : get_IPQTOQ,get_NZ

  implicit none

  integer, intent(in) :: mat_dim
  real(8), intent(in) :: dmat_0(*)
  real(8), intent(in) :: fmat(*)
! fmat should be inout

  !miro: control print out of fmat at the beginnig
  !print *,'xcint_potential_rks: entering fmat matrix, get_nz()=',get_nz()
  !call prqmat(fmat,mat_dim,mat_dim,mat_dim,mat_dim, get_nz(), get_ipqtoq(1,0), 6)

  call integrate_xc(xc_mat_dim      = mat_dim, &
                    xc_nz           = 1,       &
                    xc_dmat_0       = dmat_0,  &
                    xc_nr_dmat      = 0,       &
                    xc_nr_fmat      = 1,       &
                    xc_fmat         = fmat,    &
                    xc_do_potential = .true.)

 ! print *,'xcint_potential_rks: exit fmat matrix'
 ! call prqmat(fmat,mat_dim,mat_dim,mat_dim,mat_dim, get_nz(), get_ipqtoq(1,0), 6)

end subroutine

subroutine xcint_launch_slave_process()

  use xcint_main

  implicit none

  call integrate_xc(xc_mat_dim = 1,         &
                    xc_nz      = 1,         &
                    xc_dmat_0  = (/0.0d0/), &
                    xc_nr_dmat = 0,         &
                    xc_nr_fmat = 0)

end subroutine

function xcint_get_xc_energy()
  use xcint_main
  implicit none
  real(8) :: xcint_get_xc_energy
  xcint_get_xc_energy = get_xc_energy()
end function

function xcint_get_xc_mat_energy()
  use xcint_main
  implicit none
  real(8) :: xcint_get_xc_mat_energy
  xcint_get_xc_mat_energy = get_xc_mat_energy()
end function

function xcint_get_nr_electrons_integrated()
  use xcint_main
  implicit none
  real(8) :: xcint_get_nr_electrons_integrated
  xcint_get_nr_electrons_integrated = get_nr_electrons_integrated()
end function

subroutine xcint_parse_functional(line)

  use xc_derv
  use interface_functional_read

  implicit none

  character(80), intent(in) :: line
  real(8)                   :: hfx_out, mu_out, beta_out

  call parse_functional(line, xc_fun, hfx_out, mu_out, beta_out, .true.)
  fun_is_automatic = .true.

end subroutine
