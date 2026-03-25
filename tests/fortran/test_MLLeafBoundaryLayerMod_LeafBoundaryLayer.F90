program test_MLLeafBoundaryLayerMod_LeafBoundaryLayer
  !
  ! Test driver for MLLeafBoundaryLayerMod :: CalcLeafBoundaryLayer.
  !
  ! Reads scalar inputs from a namelist on stdin, calls CalcLeafBoundaryLayer,
  ! and writes KEY=VALUE pairs to stdout for parsing by Python.
  !
  use shr_kind_mod,           only : r8 => shr_kind_r8
  use clm_varctl,             only : iulog
  use MLLeafBoundaryLayerMod, only : CalcLeafBoundaryLayer
  implicit none

  real(r8) :: d           ! Leaf dimension (m)
  real(r8) :: u           ! Wind speed (m/s)
  real(r8) :: tleaf       ! Leaf temperature (K)
  real(r8) :: tair        ! Air temperature at leaf level (K)
  real(r8) :: tref        ! Air temperature at reference height (K)
  real(r8) :: pref        ! Air pressure at reference height (Pa)
  real(r8) :: rhomol      ! Molar density of air (mol/m3)
  integer  :: gb_type_in  ! Convection regime (1/2/3)

  real(r8) :: gbh, gbv, gbc

  namelist /inputs/ d, u, tleaf, tair, tref, pref, rhomol, gb_type_in

  iulog = 6

  ! Defaults
  d          = 0.04_r8
  u          = 2.0_r8
  tleaf      = 300.0_r8
  tair       = 298.0_r8
  tref       = 298.0_r8
  pref       = 101325.0_r8
  rhomol     = 41.5_r8
  gb_type_in = 2

  read(*, nml=inputs)

  call CalcLeafBoundaryLayer(d, u, tleaf, tair, tref, pref, rhomol, gb_type_in, gbh, gbv, gbc)

  write(*, '(A,ES24.16)') 'gbh=', gbh
  write(*, '(A,ES24.16)') 'gbv=', gbv
  write(*, '(A,ES24.16)') 'gbc=', gbc

end program test_MLLeafBoundaryLayerMod_LeafBoundaryLayer
