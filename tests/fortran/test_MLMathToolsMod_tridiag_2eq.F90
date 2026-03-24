! Test driver for MLMathToolsMod :: tridiag_2eq
!
! Solves two coupled tridiagonal systems for air temperature (t) and
! water vapor (q) at each canopy layer:
!
!   a1(i)*t(i-1) + b11(i)*t(i) + b12(i)*q(i) + c1(i)*t(i+1) = d1(i)
!   a2(i)*q(i-1) + b21(i)*t(i) + b22(i)*q(i) + c2(i)*q(i+1) = d2(i)
!
! Note: a*(1) and c*(n) are not referenced.
!
! Input namelist (stdin):
!   &inputs
!     n                        = <integer>   (system size, max 20)
!     a1,b11,b12,c1,d1         = <n floats each>
!     a2,b21,b22,c2,d2         = <n floats each>
!   /
!
! Output (stdout):
!   t_1=<value>
!   t_2=<value>
!   ...
!   q_1=<value>
!   q_2=<value>
!   ...

program test_tridiag_2eq

  use shr_kind_mod, only : r8 => shr_kind_r8
  use clm_varctl,   only : iulog
  use MLclm_varpar, only : nlevmlcan
  use MLMathToolsMod, only : tridiag_2eq

  implicit none

  integer :: n, i

  ! Declare at full module size to satisfy bounds-checking compiler flags
  real(r8) :: a1(nlevmlcan), b11(nlevmlcan), b12(nlevmlcan), c1(nlevmlcan), d1(nlevmlcan)
  real(r8) :: a2(nlevmlcan), b21(nlevmlcan), b22(nlevmlcan), c2(nlevmlcan), d2(nlevmlcan)
  real(r8) :: t(nlevmlcan), q(nlevmlcan)

  namelist /inputs/ n, a1, b11, b12, c1, d1, a2, b21, b22, c2, d2

  iulog = 6

  ! Defaults: identity system for both equations (no coupling)
  n   = 1
  a1  = 0._r8; b11 = 1._r8; b12 = 0._r8; c1  = 0._r8; d1  = 0._r8
  a2  = 0._r8; b21 = 0._r8; b22 = 1._r8; c2  = 0._r8; d2  = 0._r8

  read(*, nml=inputs)

  call tridiag_2eq(a1, b11, b12, c1, d1, a2, b21, b22, c2, d2, t, q, n)

  do i = 1, n
    write(*, '(A,I0,A,ES24.16)') 't_', i, '=', t(i)
    write(*, '(A,I0,A,ES24.16)') 'q_', i, '=', q(i)
  end do

end program test_tridiag_2eq
