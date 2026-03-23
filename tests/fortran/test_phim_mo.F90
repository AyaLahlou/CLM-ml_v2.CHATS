! Test driver for MLCanopyTurbulenceMod :: phim_monin_obukhov and phic_monin_obukhov
!
! Monin-Obukhov phi stability correction functions (Bonan et al. 2018, eq. A10-A11):
!
!   phim(zeta):
!     stable   (zeta >= 0): 1 + 5*zeta
!     unstable (zeta <  0): (1 - 16*zeta)^(-1/4)
!
!   phic(zeta):
!     stable   (zeta >= 0): 1 + 5*zeta
!     unstable (zeta <  0): (1 - 16*zeta)^(-1/2)
!
! Reference values:
!   zeta=0:  phim=1, phic=1  (neutral)
!   zeta=1:  phim=6, phic=6  (stable)
!
! Input namelist (stdin):
!   &inputs  zeta=<real>  /
!
! Output (stdout):
!   phim=<value>
!   phic=<value>

program test_phim_mo

  use shr_kind_mod, only : r8 => shr_kind_r8
  use clm_varctl, only : iulog
  use MLCanopyTurbulenceMod, only : phim_monin_obukhov, phic_monin_obukhov

  implicit none

  real(r8) :: zeta, phim, phic
  namelist /inputs/ zeta

  iulog = 6

  zeta = 0._r8
  read(*, nml=inputs)

  phim = phim_monin_obukhov(zeta)
  phic = phic_monin_obukhov(zeta)

  write(*, '(A,ES24.16)') 'phim=', phim
  write(*, '(A,ES24.16)') 'phic=', phic

end program test_phim_mo
