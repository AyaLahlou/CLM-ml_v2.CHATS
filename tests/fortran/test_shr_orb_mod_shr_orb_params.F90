! Test driver for shr_orb_mod :: shr_orb_params
!
! Computes orbital parameters for a given calendar year.
!
! Input namelist (stdin):
!   &inputs  iyear_AD=<int>  /
!     iyear_AD - calendar year (e.g. 2000 for year 2000 AD)
!
! Output (stdout):
!   eccen=<value>   orbital eccentricity (dimensionless, ~0.0167 today)
!   obliq=<value>   obliquity in degrees (~23.44 today)
!   mvelp=<value>   moving vernal equinox longitude of perihelion (degrees)
!   obliqr=<value>  obliquity in radians
!   lambm0=<value>  mean longitude of perihelion at vernal equinox (radians)
!   mvelpp=<value>  mvelp in radians plus pi

program test_shr_orb_params

  use shr_kind_mod, only : r8 => shr_kind_r8
  use clm_varctl, only : iulog
  use shr_orb_mod, only : shr_orb_params

  implicit none

  integer  :: iyear_AD
  real(r8) :: eccen, obliq, mvelp, obliqr, lambm0, mvelpp
  namelist /inputs/ iyear_AD

  iulog = 6

  iyear_AD = 2000
  read(*, nml=inputs)

  call shr_orb_params(iyear_AD, eccen, obliq, mvelp, obliqr, lambm0, mvelpp)

  write(*, '(A,ES24.16)') 'eccen=',  eccen
  write(*, '(A,ES24.16)') 'obliq=',  obliq
  write(*, '(A,ES24.16)') 'mvelp=',  mvelp
  write(*, '(A,ES24.16)') 'obliqr=', obliqr
  write(*, '(A,ES24.16)') 'lambm0=', lambm0
  write(*, '(A,ES24.16)') 'mvelpp=', mvelpp

end program test_shr_orb_params
