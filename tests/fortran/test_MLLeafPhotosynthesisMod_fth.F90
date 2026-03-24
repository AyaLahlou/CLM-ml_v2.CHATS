! Test driver for MLLeafPhotosynthesisMod :: fth
!
! High-temperature inhibition of photosynthesis:
!   fth(tl, hd, se, c) = c / (1 + exp((-hd + se*tl) / (rgas*tl)))
!
! When called with c = fth25(hd, se), this is normalised so that fth(298.15) = 1.
!
! Input namelist (stdin):
!   &inputs  tl=<real>  hd=<real>  se=<real>  c=<real>  /
!     tl  - leaf temperature (K)
!     hd  - deactivation energy (J/mol), e.g. 150000
!     se  - entropy term (J/mol/K), e.g. 490
!     c   - scaling factor; set to fth25(hd,se) for normalised response
!
! Output (stdout):
!   ans=<value>

program test_fth

  use shr_kind_mod, only : r8 => shr_kind_r8
  use clm_varctl, only : iulog
  use MLLeafPhotosynthesisMod, only : fth

  implicit none

  real(r8) :: tl, hd, se, c, ans
  namelist /inputs/ tl, hd, se, c

  iulog = 6

  tl = 298.15_r8; hd = 150000._r8; se = 490._r8; c = 1._r8
  read(*, nml=inputs)

  ans = fth(tl, hd, se, c)

  write(*, '(A,ES24.16)') 'ans=', ans

end program test_fth
