! Test driver for MLMathToolsMod :: beta_function
!
! Returns B(a, b) = Gamma(a)*Gamma(b) / Gamma(a+b).
!
! Input namelist (stdin):
!   &inputs  a=<real>  b=<real>  /
!
! Output (stdout):
!   beta=<value>

program test_beta_function

  use shr_kind_mod, only : r8 => shr_kind_r8
  use clm_varctl, only : iulog
  use MLMathToolsMod, only : beta_function

  implicit none

  real(r8) :: a, b, beta
  namelist /inputs/ a, b

  iulog = 6

  a = 1._r8; b = 1._r8
  read(*, nml=inputs)

  beta = beta_function(a, b)

  write(*, '(A,ES24.16)') 'beta=', beta

end program test_beta_function
