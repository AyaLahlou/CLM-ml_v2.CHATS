! Test driver for MLCanopyTurbulenceMod :: GetPsiRSL
!
! Computes the RSL-modified stability functions psi for momentum and
! scalars between height za and the canopy top hc.
! See Bonan et al. (2018), appendix A2.
!
! Formula:
!   psim = -psim(za) + psim(hc) + psihat_m(za) - psihat_m(hc) + vkc/beta
!   psic = -psic(za) + psic(hc) + psihat_c(za) - psihat_c(hc)
!
! Key property at za = hc (both Monin-Obukhov and RSL terms cancel):
!   psim = vkc / beta   (= 0.4 / beta)
!   psic = 0.0
!
! IMPORTANT: LookupPsihatINI must be called before GetPsiRSL to initialise
!            the RSL psihat look-up table grids.
!
! Input namelist (stdin):
!   &inputs  za=<real>  hc=<real>  disp=<real>  obu=<real>
!            beta=<real>  PrSc=<real>  /
!
! Output (stdout):
!   psim=<value>        psi for momentum including RSL influence
!   psic=<value>        psi for scalars including RSL influence
!   psim2=<value>       Monin-Obukhov psi_m evaluated at hc
!   psim_hat2=<value>   RSL psihat_m evaluated at hc

program test_GetPsiRSL

  use shr_kind_mod,          only : r8 => shr_kind_r8
  use clm_varctl,            only : iulog
  use MLCanopyTurbulenceMod, only : GetPsiRSL, LookupPsihatINI

  implicit none

  real(r8) :: za, hc, disp, obu, beta, PrSc
  real(r8) :: psim, psic, psim2, psim_hat2

  namelist /inputs/ za, hc, disp, obu, beta, PrSc

  iulog = 6

  ! Defaults: reference height 30 m above 20 m canopy, moderately unstable
  za   = 30._r8
  hc   = 20._r8
  disp = 12._r8
  obu  = -100._r8
  beta = 0.3_r8
  PrSc = 0.5_r8

  read(*, nml=inputs)

  ! Initialise RSL psihat look-up tables (required before calling GetPsiRSL)
  call LookupPsihatINI()

  call GetPsiRSL(za, hc, disp, obu, beta, PrSc, psim, psic, psim2, psim_hat2)

  write(*, '(A,ES24.16)') 'psim=',      psim
  write(*, '(A,ES24.16)') 'psic=',      psic
  write(*, '(A,ES24.16)') 'psim2=',     psim2
  write(*, '(A,ES24.16)') 'psim_hat2=', psim_hat2

end program test_GetPsiRSL
