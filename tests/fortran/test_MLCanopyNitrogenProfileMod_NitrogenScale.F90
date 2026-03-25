program test_MLCanopyNitrogenProfileMod_NitrogenScale
  !
  ! Test driver for MLCanopyNitrogenProfileMod :: CalcNitrogenScale.
  !
  ! Reads scalar inputs from a namelist on stdin, calls CalcNitrogenScale,
  ! and writes KEY=VALUE pairs to stdout for parsing by Python.
  !
  use shr_kind_mod,                only : r8 => shr_kind_r8
  use clm_varctl,                  only : iulog
  use MLCanopyNitrogenProfileMod,  only : CalcNitrogenScale
  implicit none

  real(r8) :: kn                  ! Leaf nitrogen decay coefficient
  real(r8) :: pai_above           ! Cumulative PAI above this layer
  real(r8) :: dpai                ! Layer plant area index (m2/m2)
  real(r8) :: kb                  ! Direct beam extinction coefficient
  real(r8) :: clump_fac           ! Foliage clumping index (-)
  real(r8) :: fracsun             ! Sunlit fraction of layer (-)
  real(r8) :: tbi                 ! Cumulative direct beam transmittance (-)
  integer  :: leaf_optics_type_in ! 0 = Bonan (2021); 1 = thin-layer

  real(r8) :: nscale_sun, nscale_sha

  namelist /inputs/ kn, pai_above, dpai, kb, clump_fac, fracsun, tbi, leaf_optics_type_in

  iulog = 6

  ! Defaults
  kn                  = 0.3_r8
  pai_above           = 0.0_r8
  dpai                = 0.5_r8
  kb                  = 0.5_r8
  clump_fac           = 0.8_r8
  fracsun             = 0.5_r8
  tbi                 = 1.0_r8
  leaf_optics_type_in = 0

  read(*, nml=inputs)

  call CalcNitrogenScale(kn, pai_above, dpai, kb, clump_fac, fracsun, tbi, &
                         leaf_optics_type_in, nscale_sun, nscale_sha)

  write(*, '(A,ES24.16)') 'nscale_sun=', nscale_sun
  write(*, '(A,ES24.16)') 'nscale_sha=', nscale_sha

end program test_MLCanopyNitrogenProfileMod_NitrogenScale
