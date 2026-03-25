! Test driver for MLCanopyWaterMod :: CalcWettedFraction
!
! Computes wetted (fwet) and dry (fdry) fractions for a single canopy layer.
!
! Formulas (using module constants from MLclm_varcon):
!   h2ocanmx = dewmx * dpai           (dewmx = 0.1 kg H2O/m2 leaf)
!   fwet = clamp((h2ocan/h2ocanmx)^fwet_exponent, 0, max_fwet)
!                                     (fwet_exponent = 0.67,
!                                      max_fwet = maximum_leaf_wetted_fraction = 0.05)
!   fdry = (1 - fwet) * dlai/dpai
!
! Special cases:
!   dpai <= 0         : fwet = fdry = 0
!   h2ocan = 0        : fwet = 0,   fdry = dlai/dpai
!   h2ocan >= h2ocanmx: fwet hits the cap (0.05)
!
! Input namelist (stdin):
!   &inputs  h2ocan=<real>  dpai=<real>  dlai=<real>  /
!
! Output (stdout):
!   fwet=<value>
!   fdry=<value>

program test_WettedFraction

  use shr_kind_mod,      only : r8 => shr_kind_r8
  use clm_varctl,        only : iulog
  use MLCanopyWaterMod,  only : CalcWettedFraction

  implicit none

  real(r8) :: h2ocan, dpai, dlai
  real(r8) :: fwet, fdry

  namelist /inputs/ h2ocan, dpai, dlai

  iulog = 6

  h2ocan = 0._r8
  dpai   = 1._r8
  dlai   = 0.8_r8

  read(*, nml=inputs)

  call CalcWettedFraction(h2ocan, dpai, dlai, fwet, fdry)

  write(*, '(A,ES24.16)') 'fwet=', fwet
  write(*, '(A,ES24.16)') 'fdry=', fdry

end program test_WettedFraction
