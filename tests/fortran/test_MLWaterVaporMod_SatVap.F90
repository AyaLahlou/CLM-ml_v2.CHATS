! Test driver for MLWaterVaporMod :: SatVap
!
! Computes saturation vapor pressure (es, Pa) and its temperature
! derivative (desdt, Pa/K) using the polynomial fit of Flatau et al. (1992).
!
! Valid temperature range: -75°C to +100°C (clamped outside this range).
! Water branch used for T >= 273.15 K; ice branch for T < 273.15 K.
!
! Input namelist (stdin):
!   &inputs  t=<real>  /   (temperature in K)
!
! Output (stdout):
!   es=<value>      (Pa)
!   desdt=<value>   (Pa/K)

program test_SatVap

  use shr_kind_mod, only : r8 => shr_kind_r8
  use clm_varctl, only : iulog
  use MLWaterVaporMod, only : SatVap

  implicit none

  real(r8) :: t, es, desdt
  namelist /inputs/ t

  iulog = 6

  t = 273.15_r8
  read(*, nml=inputs)

  call SatVap(t, es, desdt)

  write(*, '(A,ES24.16)') 'es=', es
  write(*, '(A,ES24.16)') 'desdt=', desdt

end program test_SatVap
