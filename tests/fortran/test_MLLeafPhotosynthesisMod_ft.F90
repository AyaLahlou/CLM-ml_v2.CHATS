! Test driver for MLLeafPhotosynthesisMod :: ft
!
! Arrhenius-type photosynthesis temperature response, normalised to 1.0 at 25°C:
!   ft(tl, ha) = exp( ha / (rgas * 298.15) * (1 - 298.15 / tl) )
!
! Identity: ft(298.15, ha) == 1.0 for any ha.
!
! Input namelist (stdin):
!   &inputs  tl=<real>  ha=<real>  /
!     tl  - leaf temperature (K)
!     ha  - activation energy (J/mol), e.g. 65330 (vcmax), 43540 (jmax)
!
! Output (stdout):
!   ans=<value>

program test_ft

  use shr_kind_mod, only : r8 => shr_kind_r8
  use clm_varctl, only : iulog
  use MLLeafPhotosynthesisMod, only : ft

  implicit none

  real(r8) :: tl, ha, ans
  namelist /inputs/ tl, ha

  iulog = 6

  tl = 298.15_r8; ha = 65330._r8
  read(*, nml=inputs)

  ans = ft(tl, ha)

  write(*, '(A,ES24.16)') 'ans=', ans

end program test_ft
