! Test driver for MLMathToolsMod :: beta_distribution_pdf
!
! Returns the beta distribution PDF: f(x; a, b).
!
! Input namelist (stdin):
!   &inputs  a=<real>  b=<real>  x=<real>  /
!
! Output (stdout):
!   beta_pdf=<value>

program test_beta_pdf

  use shr_kind_mod, only : r8 => shr_kind_r8
  use clm_varctl, only : iulog
  use MLMathToolsMod, only : beta_distribution_pdf

  implicit none

  real(r8) :: a, b, x, beta_pdf
  namelist /inputs/ a, b, x

  iulog = 6

  a = 2._r8; b = 2._r8; x = 0.5_r8
  read(*, nml=inputs)

  beta_pdf = beta_distribution_pdf(a, b, x)

  write(*, '(A,ES24.16)') 'beta_pdf=', beta_pdf

end program test_beta_pdf
