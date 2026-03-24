! Test driver for MLWaterVaporMod :: LatVap
!
! Returns the molar latent heat of vaporisation (J/mol).
! Uses hvap (evaporation) for T > 273.15 K, hsub (sublimation) otherwise.
! Both are multiplied by the molar mass of water (mmh2o = 0.01802 kg/mol).
!
! Input namelist (stdin):
!   &inputs  t=<real>  /   (temperature in K)
!
! Output (stdout):
!   lambda=<value>   (J/mol)

program test_LatVap

  use shr_kind_mod, only : r8 => shr_kind_r8
  use clm_varctl, only : iulog
  use MLWaterVaporMod, only : LatVap

  implicit none

  real(r8) :: t, lambda
  namelist /inputs/ t

  iulog = 6

  t = 273.15_r8
  read(*, nml=inputs)

  lambda = LatVap(t)

  write(*, '(A,ES24.16)') 'lambda=', lambda

end program test_LatVap
