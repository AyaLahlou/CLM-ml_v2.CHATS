! Test driver for MLCanopyTurbulenceMod :: psim_monin_obukhov and psic_monin_obukhov
!
! Monin-Obukhov psi integrated stability correction functions (Bonan et al. 2018, eq. A12-A13):
!
!   psim(zeta):
!     stable   (zeta >= 0): -5*zeta
!     unstable (zeta <  0): 2*ln((1+x)/2) + ln((1+x^2)/2) - 2*atan(x) + pi/2
!                           where x = (1-16*zeta)^(1/4)
!
!   psic(zeta):
!     stable   (zeta >= 0): -5*zeta
!     unstable (zeta <  0): 2*ln((1+x^2)/2)
!                           where x = (1-16*zeta)^(1/4)
!
! Reference values:
!   zeta=0:  psim=0, psic=0  (neutral, no correction)
!   zeta=1:  psim=-5, psic=-5  (stable)
!
! Input namelist (stdin):
!   &inputs  zeta=<real>  /
!
! Output (stdout):
!   psim=<value>
!   psic=<value>

program test_psim_mo

  use shr_kind_mod, only : r8 => shr_kind_r8
  use clm_varctl, only : iulog
  use MLCanopyTurbulenceMod, only : psim_monin_obukhov, psic_monin_obukhov

  implicit none

  real(r8) :: zeta, psim, psic
  namelist /inputs/ zeta

  iulog = 6

  zeta = 0._r8
  read(*, nml=inputs)

  psim = psim_monin_obukhov(zeta)
  psic = psic_monin_obukhov(zeta)

  write(*, '(A,ES24.16)') 'psim=', psim
  write(*, '(A,ES24.16)') 'psic=', psic

end program test_psim_mo
