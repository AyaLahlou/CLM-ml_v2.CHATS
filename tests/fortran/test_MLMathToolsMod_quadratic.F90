! Test driver for MLMathToolsMod :: quadratic
!
! Solves ax^2 + bx + c = 0 for its two roots r1, r2.
!
! Input namelist (stdin):
!   &inputs  a=<real>  b=<real>  c=<real>  /
!
! Output (stdout):
!   r1=<value>
!   r2=<value>

program test_quadratic

  use shr_kind_mod, only : r8 => shr_kind_r8
  use clm_varctl, only : iulog
  use MLMathToolsMod, only : quadratic

  implicit none

  real(r8) :: a, b, c, r1, r2
  namelist /inputs/ a, b, c

  ! Route Fortran log output to stdout
  iulog = 6

  a = 0._r8; b = 0._r8; c = 0._r8
  read(*, nml=inputs)

  call quadratic(a, b, c, r1, r2)

  write(*, '(A,ES24.16)') 'r1=', r1
  write(*, '(A,ES24.16)') 'r2=', r2

end program test_quadratic
