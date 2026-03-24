! Test driver for MLCanopyTurbulenceMod :: GetPrSc
!
! Calculates the effective Prandtl number (Pr) / Schmidt number (Sc) at the
! canopy top.  The formula is (Bonan et al. 2018, eqs. A25, A34):
!
!   PrSc = Pr0 + Pr1 * tanh(Pr2 * LcL)
!
! where the constants from MLclm_varcon default to:
!   Pr0 = 0.5   (neutral value)
!   Pr1 = 0.3   (amplitude of stability variation)
!   Pr2 = 2.0   (scale of stability variation)
!
! When sparse_canopy_type = 1 the result is blended toward 1.0:
!   PrSc = (1 - beta_neutral/beta_neutral_max) + (beta_neutral/beta_neutral_max)*PrSc
!
! Input namelist (stdin):
!   &inputs  beta_neutral=<real>  LcL=<real>  sparse_canopy_type_in=<int>  /
!
! Output (stdout):
!   PrSc=<value>

program test_GetPrSc

  use shr_kind_mod, only : r8 => shr_kind_r8
  use clm_varctl,   only : iulog
  use MLclm_varctl, only : sparse_canopy_type
  use MLclm_varcon, only : beta_neutral_max
  use MLCanopyTurbulenceMod, only : GetPrSc

  implicit none

  real(r8) :: beta_neutral, LcL, PrSc
  integer  :: sparse_canopy_type_in

  namelist /inputs/ beta_neutral, LcL, sparse_canopy_type_in

  iulog = 6

  beta_neutral          = 0.3_r8
  LcL                   = 0._r8
  sparse_canopy_type_in = 0

  read(*, nml=inputs)
  sparse_canopy_type = sparse_canopy_type_in

  call GetPrSc(beta_neutral, beta_neutral_max, LcL, PrSc)

  write(*, '(A,ES24.16)') 'PrSc=', PrSc

end program test_GetPrSc
