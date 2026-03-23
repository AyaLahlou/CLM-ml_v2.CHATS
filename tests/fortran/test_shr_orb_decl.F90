! Test driver for shr_orb_mod :: shr_orb_decl
!
! Computes solar declination (delta) and Earth-sun distance factor (eccf)
! for a given calendar day, given orbital parameters.
!
! Workflow: call shr_orb_params to get orbital parameters for the year,
! then call shr_orb_decl to get declination for a specific day.
!
! Input namelist (stdin):
!   &inputs  iyear_AD=<int>  calday=<real>  /
!     iyear_AD - calendar year (e.g. 2007)
!     calday   - day of year including fraction (1.0 = 1 Jan, 80.0 ~ vernal equinox)
!
! Output (stdout):
!   eccen=<value>    orbital eccentricity
!   obliqr=<value>   obliquity in radians
!   delta=<value>    solar declination (radians); ~0 at equinoxes, +obliqr at summer solstice
!   eccf=<value>     Earth-sun distance factor (1/r)^2; >1 when closer to sun

program test_shr_orb_decl

  use shr_kind_mod, only : r8 => shr_kind_r8
  use clm_varctl, only : iulog
  use shr_orb_mod, only : shr_orb_params, shr_orb_decl

  implicit none

  integer  :: iyear_AD
  real(r8) :: calday
  real(r8) :: eccen, obliq, mvelp, obliqr, lambm0, mvelpp
  real(r8) :: delta, eccf
  namelist /inputs/ iyear_AD, calday

  iulog = 6

  iyear_AD = 2007; calday = 80._r8
  read(*, nml=inputs)

  ! Get orbital parameters for this year
  call shr_orb_params(iyear_AD, eccen, obliq, mvelp, obliqr, lambm0, mvelpp)

  ! Get solar declination for this calendar day
  call shr_orb_decl(calday, eccen, mvelpp, lambm0, obliqr, delta, eccf)

  write(*, '(A,ES24.16)') 'eccen=',  eccen
  write(*, '(A,ES24.16)') 'obliqr=', obliqr
  write(*, '(A,ES24.16)') 'delta=',  delta
  write(*, '(A,ES24.16)') 'eccf=',   eccf

end program test_shr_orb_decl
