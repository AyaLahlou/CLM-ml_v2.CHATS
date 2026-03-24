! Test driver for MLCanopyTurbulenceMod :: GetBeta
!
! Calculates beta = u*/u(h), the ratio of friction velocity to wind speed
! at the canopy top, for a given stability state.
! See Bonan et al. (2018), eqs. (A22)-(A24).
!
! Stability is parameterised by:
!   LcL = Lc / obu   (canopy density length scale / Obukhov length)
!     LcL = 0  → neutral conditions
!     LcL < 0  → unstable
!     LcL > 0  → stable
!
! Defining property: beta * phim(LcL * beta^2) = beta_neutral
!
! Input namelist (stdin):
!   &inputs  beta_neutral=<real>  LcL=<real>  /
!
! Output (stdout):
!   beta=<value>

program test_GetBeta

  use shr_kind_mod, only : r8 => shr_kind_r8
  use clm_varctl,   only : iulog
  use MLCanopyTurbulenceMod, only : GetBeta

  implicit none

  real(r8) :: beta_neutral, LcL, beta

  namelist /inputs/ beta_neutral, LcL

  iulog = 6

  beta_neutral = 0.3_r8
  LcL          = 0._r8

  read(*, nml=inputs)

  call GetBeta(beta_neutral, LcL, beta)

  write(*, '(A,ES24.16)') 'beta=', beta

end program test_GetBeta
