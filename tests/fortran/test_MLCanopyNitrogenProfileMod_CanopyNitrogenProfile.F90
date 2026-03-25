program test_MLCanopyNitrogenProfileMod_CanopyNitrogenProfile
  !
  ! Parent-level integration test for the CanopyNitrogenProfile loop.
  !
  ! Replicates the CanopyNitrogenProfile multi-layer loop using the public
  ! CalcNitrogenScale helper, then validates the key consistency condition:
  !
  !   sum_ic(vcmax25_profile(ic) * dpai(ic)) = vcmax25top * (1 - exp(-kn*totalPAI)) / kn
  !
  ! This is the same assertion that CanopyNitrogenProfile makes internally.
  ! The test exercises the full Beer's-law nitrogen integration across all canopy
  ! layers, not just a single layer as test_NitrogenScale.exe does.
  !
  ! Conventions (matching CanopyNitrogenProfile):
  !   ic = ncan  →  top-most layer  (pai_above = 0 initially)
  !   ic = 1     →  bottom layer    (largest pai_above)
  !   totalPAI   = lai + sai
  !
  use shr_kind_mod,               only : r8 => shr_kind_r8
  use clm_varctl,                 only : iulog
  use MLCanopyNitrogenProfileMod, only : CalcNitrogenScale
  implicit none

  integer, parameter :: ncan_max = 10

  integer  :: ncan, leaf_optics_type_in, ic
  real(r8) :: dpai(ncan_max)    ! Layer plant area index (m2/m2); ic=ncan is top
  real(r8) :: fracsun(ncan_max) ! Sunlit fraction per layer (-)
  real(r8) :: kb(ncan_max)      ! Direct beam extinction coefficient per layer (-)
  real(r8) :: tbi(ncan_max)     ! Cumulative direct beam transmittance onto each layer (-)
  real(r8) :: clump_fac         ! Foliage clumping index (-)
  real(r8) :: kn                ! Leaf nitrogen decay coefficient (-)
  real(r8) :: vcmax25top        ! Canopy top Vcmax at 25 C (umol/m2/s)
  real(r8) :: lai               ! Leaf area index (m2/m2)
  real(r8) :: sai               ! Stem area index (m2/m2)

  real(r8) :: nscale_sun(ncan_max), nscale_sha(ncan_max), vcmax25_profile(ncan_max)
  real(r8) :: pai_above, numerical, analytical

  namelist /inputs/ ncan, dpai, fracsun, kb, tbi, clump_fac, kn, vcmax25top, lai, sai, &
                    leaf_optics_type_in

  iulog = 6

  ! --- Defaults: 3-layer uniform canopy with kn=0.3 ---
  ! Convention: index 3 = top layer, index 1 = bottom layer
  ncan                = 3
  dpai(:)             = 0._r8
  dpai(1:3)           = 0.5_r8
  fracsun(:)          = 0._r8
  fracsun(1)          = 0.20_r8      ! bottom
  fracsun(2)          = 0.35_r8      ! middle
  fracsun(3)          = 0.50_r8      ! top
  kb(:)               = 0._r8
  kb(1:3)             = 0.5_r8
  tbi(:)              = 0._r8
  tbi(3)              = 1.0_r8       ! top layer: no attenuation above
  tbi(2)              = 0.7788007831_r8  ! exp(-0.5*1.0*0.5)
  tbi(1)              = 0.6065306597_r8  ! exp(-0.5*1.0*1.0)
  clump_fac           = 1.0_r8
  kn                  = 0.3_r8
  vcmax25top          = 60.0_r8
  lai                 = 1.5_r8
  sai                 = 0.0_r8
  leaf_optics_type_in = 0

  read(*, nml=inputs)

  ! Initialise outputs
  nscale_sun(:)       = 0._r8
  nscale_sha(:)       = 0._r8
  vcmax25_profile(:)  = 0._r8

  ! Replicate the CanopyNitrogenProfile loop (top to bottom)
  pai_above = 0._r8
  do ic = ncan, 1, -1
    if (dpai(ic) > 0._r8) then
      call CalcNitrogenScale(kn, pai_above, dpai(ic), kb(ic), clump_fac, &
                             fracsun(ic), tbi(ic), leaf_optics_type_in,   &
                             nscale_sun(ic), nscale_sha(ic))
      vcmax25_profile(ic) = vcmax25top * ( nscale_sun(ic) * fracsun(ic)         &
                                         + nscale_sha(ic) * (1._r8 - fracsun(ic)) )
    end if
    pai_above = pai_above + dpai(ic)
  end do

  ! Integration sums – the main consistency check
  numerical  = sum(vcmax25_profile(1:ncan) * dpai(1:ncan))
  analytical = vcmax25top * (1._r8 - exp(-kn * (lai + sai))) / kn

  ! Scalar summary outputs
  write(*, '(A,ES24.16)') 'numerical=',  numerical
  write(*, '(A,ES24.16)') 'analytical=', analytical

  ! Per-layer outputs (1 = bottom, ncan = top)
  do ic = 1, ncan
    write(*, '(A,I2.2,A,ES24.16)') 'nscale_sun_',      ic, '=', nscale_sun(ic)
    write(*, '(A,I2.2,A,ES24.16)') 'nscale_sha_',      ic, '=', nscale_sha(ic)
    write(*, '(A,I2.2,A,ES24.16)') 'vcmax25_profile_', ic, '=', vcmax25_profile(ic)
  end do

end program test_MLCanopyNitrogenProfileMod_CanopyNitrogenProfile
