! Test driver for MLMathToolsMod :: log_gamma_function
!
! Returns ln(Gamma(x)).
!
! Input namelist (stdin):
!   &inputs  x=<real>  /
!
! Output (stdout):
!   gammaln=<value>

program test_log_gamma

  use shr_kind_mod, only : r8 => shr_kind_r8
  use clm_varctl, only : iulog
  use MLMathToolsMod, only : log_gamma_function

  implicit none

  real(r8) :: x, gammaln
  namelist /inputs/ x

  iulog = 6

  x = 1._r8
  read(*, nml=inputs)

  gammaln = log_gamma_function(x)

  write(*, '(A,ES24.16)') 'gammaln=', gammaln

end program test_log_gamma
