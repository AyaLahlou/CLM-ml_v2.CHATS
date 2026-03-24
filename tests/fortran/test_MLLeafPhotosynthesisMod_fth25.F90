! Test driver for MLLeafPhotosynthesisMod :: fth25
!
! Scaling factor for photosynthesis temperature inhibition at 25°C (298.15 K):
!   fth25(hd, se) = 1 + exp((-hd + se * 298.15) / (rgas * 298.15))
!
! This value is always > 1 for typical parameter values where hd > se * 298.15.
! It is passed as argument c to fth() so that fth(298.15) = 1.
!
! Input namelist (stdin):
!   &inputs  hd=<real>  se=<real>  /
!     hd  - deactivation energy (J/mol), e.g. 150000
!     se  - entropy term (J/mol/K), e.g. 490
!
! Output (stdout):
!   ans=<value>

program test_fth25

  use shr_kind_mod, only : r8 => shr_kind_r8
  use clm_varctl, only : iulog
  use MLLeafPhotosynthesisMod, only : fth25

  implicit none

  real(r8) :: hd, se, ans
  namelist /inputs/ hd, se

  iulog = 6

  hd = 150000._r8; se = 490._r8
  read(*, nml=inputs)

  ans = fth25(hd, se)

  write(*, '(A,ES24.16)') 'ans=', ans

end program test_fth25
