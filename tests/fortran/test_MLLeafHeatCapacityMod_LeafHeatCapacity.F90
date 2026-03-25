! Test driver for MLLeafHeatCapacityMod :: CalcLeafHeatCapacity
!
! Computes leaf heat capacity (J/K/m2 leaf) from specific leaf area (m2/gC).
! See Bonan et al. (2018) Geosci. Model Dev. eq. (A29).
!
! Formula (using module constants):
!   lma          = (1/sla) * 0.001            [kg C/m2]  (sla in m2/gC)
!   dry_weight   = lma / fcarbon              [kg DM/m2] (fcarbon = 0.5)
!   fresh_weight = dry_weight / (1 - fwater)  [kg FM/m2] (fwater = 0.7)
!   leaf_water   = fwater * fresh_weight      [kg H2O/m2]
!   cpleaf       = cpbio*dry_weight + cpliq*leaf_water   [J/K/m2]
!
! where cpbio = 4188/3 ≈ 1396 J/kg/K and cpliq = 4188 J/kg/K.
!
! The formula simplifies to:
!   cpleaf ≈ 22.336 / sla     [J/K/m2 leaf]
!
! Input namelist (stdin):
!   &inputs  sla=<real>  /    (sla in m2/gC)
!
! Output (stdout):
!   cpleaf=<value>            (J/K/m2 leaf)

program test_LeafHeatCapacity

  use shr_kind_mod,            only : r8 => shr_kind_r8
  use clm_varctl,              only : iulog
  use MLLeafHeatCapacityMod,   only : CalcLeafHeatCapacity

  implicit none

  real(r8) :: sla, cpleaf

  namelist /inputs/ sla

  iulog = 6

  sla = 0.04_r8   ! typical value: 0.04 m2/gC

  read(*, nml=inputs)

  cpleaf = CalcLeafHeatCapacity(sla)

  write(*, '(A,ES24.16)') 'cpleaf=', cpleaf

end program test_LeafHeatCapacity
