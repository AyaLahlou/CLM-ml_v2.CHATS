! Test driver for MLMathToolsMod :: beta_distribution_cdf
!
! Returns the regularised incomplete beta function I_x(a,b), i.e. the CDF
! of the beta distribution at x: F(x; a, b).
!
! Input namelist (stdin):
!   &inputs  a=<real>  b=<real>  x=<real>  /
!
! Output (stdout):
!   beta_cdf=<value>

program test_beta_cdf

  use shr_kind_mod, only : r8 => shr_kind_r8
  use clm_varctl, only : iulog
  use MLMathToolsMod, only : beta_distribution_cdf

  implicit none

  real(r8) :: a, b, x, beta_cdf
  namelist /inputs/ a, b, x

  iulog = 6

  a = 2._r8; b = 2._r8; x = 0.5_r8
  read(*, nml=inputs)

  beta_cdf = beta_distribution_cdf(a, b, x)

  write(*, '(A,ES24.16)') 'beta_cdf=', beta_cdf

end program test_beta_cdf
