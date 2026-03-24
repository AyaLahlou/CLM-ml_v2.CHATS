! Test driver for MLLeafPhotosynthesisMod :: RealizedRate
!
! Computes gross photosynthesis as the minimum or smooth co-limited rate
! from three carboxylation-limited rates:
!   ac  - Rubisco-limited rate  (umol CO2/m2 leaf/s)
!   aj  - RuBP regeneration-limited rate
!   ap  - Product-limited (C3) or CO2-limited (C4) rate
!
! Co-limitation type is controlled by colim_type in MLclm_varctl:
!   0: minimum rate    agross = min(ac, aj) for C3
!   1: smooth colimit  using quadratic with curvature parameter colim_c3a
!
! C3 vs C4 is selected by c3psn:
!   1.0 = C3 photosynthesis
!   0.0 = C4 photosynthesis
!
! Input namelist (stdin):
!   &inputs  c3psn=<real>  ac=<real>  aj=<real>  ap=<real>  colim_type_in=<int>  /
!
! Output (stdout):
!   agross=<value>   (umol CO2/m2 leaf/s)

program test_RealizedRate

  use shr_kind_mod, only : r8 => shr_kind_r8
  use clm_varctl,   only : iulog
  use MLclm_varctl, only : colim_type
  use MLLeafPhotosynthesisMod, only : RealizedRate

  implicit none

  real(r8) :: c3psn, ac, aj, ap, agross
  integer  :: colim_type_in

  namelist /inputs/ c3psn, ac, aj, ap, colim_type_in

  iulog = 6

  c3psn          = 1._r8    ! C3
  ac             = 10._r8
  aj             = 10._r8
  ap             = 10._r8
  colim_type_in  = 1

  read(*, nml=inputs)
  colim_type = colim_type_in

  call RealizedRate(c3psn, ac, aj, ap, agross)

  write(*, '(A,ES24.16)') 'agross=', agross

end program test_RealizedRate
