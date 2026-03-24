! Test driver for MLMathToolsMod :: tridiag
!
! Solves the tridiagonal system F*u = r where F is defined by vectors a, b, c.
! Matrix layout (n=3 example):
!   | b(1) c(1)   0  | | u(1) |   | r(1) |
!   | a(2) b(2) c(2) | | u(2) | = | r(2) |
!   |   0  a(3) b(3) | | u(3) |   | r(3) |
! Note: a(1) and c(n) are not referenced.
!
! Input namelist (stdin):
!   &inputs  n=<int>  a=<arr>  b=<arr>  c=<arr>  r=<arr>  /
!   (arrays must have at least n elements; max supported size is 20)
!
! Output (stdout):
!   u_1=<value>
!   u_2=<value>
!   ...

program test_tridiag

  use shr_kind_mod, only : r8 => shr_kind_r8
  use clm_varctl, only : iulog
  use MLMathToolsMod, only : tridiag

  implicit none

  integer, parameter :: nmax = 20
  integer  :: n, i
  real(r8) :: a(nmax), b(nmax), c(nmax), r(nmax), u(nmax)
  namelist /inputs/ n, a, b, c, r

  iulog = 6

  ! Defaults
  n = 3
  a = 0._r8; b = 1._r8; c = 0._r8; r = 0._r8; u = 0._r8

  read(*, nml=inputs)

  call tridiag(a, b, c, r, u, n)

  do i = 1, n
    write(*, '(A,I0,A,ES24.16)') 'u_', i, '=', u(i)
  end do

end program test_tridiag
