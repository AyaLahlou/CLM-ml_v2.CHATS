! Test driver for MLRungeKuttaMod :: RungeKuttaIni
!
! Initializes the Butcher tableau (a, b, c) for the selected Runge-Kutta method.
! runge_kutta_type and nrk are compile-time parameters in MLclm_varctl:
!   runge_kutta_type = 41  ->  4th-order Kutta method, nrk = 4
!
! Butcher tableau structure (nrk = 4, 4th-order classical Kutta):
!   c = [0, 1/2, 1/2, 1]
!   b = [1/6, 2/6, 2/6, 1/6]           (weights, sum = 1)
!   a = lower-triangular with:
!       a(2,1)=1/2
!       a(3,1)=0, a(3,2)=1/2
!       a(4,1)=0, a(4,2)=0, a(4,3)=1
!   upper triangle and diagonal: spval (unused sentinel)
!
! This driver takes no meaningful input; the namelist is accepted but ignored.
!
! Output (stdout):
!   b_1..b_nrk  =<value>
!   c_1..c_nrk  =<value>
!   a_i_j       =<value>  for i=1..nrk, j=1..nrk

program test_RungeKuttaIni

  use shr_kind_mod,   only : r8 => shr_kind_r8
  use clm_varctl,     only : iulog
  use MLclm_varctl,   only : nrk
  use MLRungeKuttaMod, only : RungeKuttaIni

  implicit none

  real(r8) :: a(nrk,nrk), b(nrk), c(nrk)
  integer  :: i, j

  ! Dummy namelist — accepts "&inputs /" from the test framework but uses no vars
  integer :: dummy = 0
  namelist /inputs/ dummy

  iulog = 6

  read(*, nml=inputs)

  call RungeKuttaIni(a, b, c)

  do i = 1, nrk
    write(*, '(A,I0,A,ES24.16)') 'b_', i, '=', b(i)
  end do
  do i = 1, nrk
    write(*, '(A,I0,A,ES24.16)') 'c_', i, '=', c(i)
  end do
  do i = 1, nrk
    do j = 1, nrk
      write(*, '(A,I0,A,I0,A,ES24.16)') 'a_', i, '_', j, '=', a(i,j)
    end do
  end do

end program test_RungeKuttaIni
