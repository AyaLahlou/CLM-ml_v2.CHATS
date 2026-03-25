module MLLeafBoundaryLayerMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Leaf boundary layer conductance
  !
  ! !USES:
  use abortutils, only : endrun
  use clm_varctl, only : iulog
  use shr_kind_mod, only : r8 => shr_kind_r8
  !
  ! !PUBLIC TYPES:
  implicit none
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: LeafBoundaryLayer
  public :: CalcLeafBoundaryLayer
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine LeafBoundaryLayer (num_filter, filter, il, mlcanopy_inst)
    !
    ! !DESCRIPTION:
    ! Leaf boundary layer conductance
    ! See Bonan (2019) Climate Change and Terrestrial Ecosystem Modeling (Chapter 10)
    !
    ! !USES:
    use clm_varcon, only : tfrz, grav
    use PatchType, only : patch
    use pftconMod, only : pftcon
    use MLclm_varcon, only : visc0, dh0, dv0, dc0, gb_factor, gbh_min
    use MLclm_varctl, only : gb_type
    use MLCanopyFluxesType, only : mlcanopy_type
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in) :: num_filter        ! Number of patches in filter
    integer, intent(in) :: filter(:)         ! Patch filter
    integer, intent(in) :: il                ! Sunlit (1) or shaded (2) leaf index
    type(mlcanopy_type), intent(inout) :: mlcanopy_inst
    !
    ! !LOCAL VARIABLES:
    integer  :: fp                           ! Filter index
    integer  :: p                            ! Patch index for CLM g/l/c/p hierarchy
    integer  :: ic                           ! Aboveground layer index
    real(r8) :: visc                         ! Kinematic viscosity (m2/s)
    real(r8) :: dh                           ! Molecular diffusivity, heat (m2/s)
    real(r8) :: dv                           ! Molecular diffusivity, H2O (m2/s)
    real(r8) :: dc                           ! Molecular diffusivity, CO2 (m2/s)
    real(r8) :: fac                          ! Correction factor for temperature and pressure
    real(r8) :: nu                           ! Nusselt number (dimensionless)
    real(r8) :: pr                           ! Prandtl number (dimensionless)
    real(r8) :: re                           ! Reynolds number (dimensionless)
    real(r8) :: gr                           ! Grashof number (dimensionless)
    real(r8) :: gbh_lam, gbv_lam, gbc_lam    ! Forced convection - laminar: conductances (mol/m2/s)
    real(r8) :: gbh_turb, gbv_turb, gbc_turb ! Forced convection - turbulent: conductances (mol/m2/s)
    real(r8) :: gbh_free, gbv_free, gbc_free ! Free convection: conductances (mol/m2/s)
    !---------------------------------------------------------------------

    associate ( &
                                                   ! *** Input ***
    dleaf     => pftcon%dleaf                 , &  ! CLM: Leaf dimension (m)
    tref      => mlcanopy_inst%tref_forcing   , &  ! Air temperature at reference height (K)
    pref      => mlcanopy_inst%pref_forcing   , &  ! Air pressure at reference height (Pa)
    rhomol    => mlcanopy_inst%rhomol_forcing , &  ! Molar density at reference height (mol/m3)
    ncan      => mlcanopy_inst%ncan_canopy    , &  ! Number of aboveground layers
    dpai      => mlcanopy_inst%dpai_profile   , &  ! Canopy layer plant area index (m2/m2)
    wind      => mlcanopy_inst%wind_profile   , &  ! Canopy layer wind speed (m/s)
    tair      => mlcanopy_inst%tair_profile   , &  ! Canopy layer air temperature (K)
    tleaf     => mlcanopy_inst%tleaf_leaf     , &  ! Leaf temperature (K)
                                                   ! *** Output ***
    gbh       => mlcanopy_inst%gbh_leaf       , &  ! Leaf boundary layer conductance: heat (mol/m2 leaf/s)
    gbv       => mlcanopy_inst%gbv_leaf       , &  ! Leaf boundary layer conductance: H2O (mol H2O/m2 leaf/s)
    gbc       => mlcanopy_inst%gbc_leaf         &  ! Leaf boundary layer conductance: CO2 (mol CO2/m2 leaf/s)
    )

    do fp = 1, num_filter
       p = filter(fp)

       ! Adjust diffusivity for temperature and pressure

       fac = 101325._r8 / pref(p) * (tref(p) / tfrz)**1.81_r8
       visc = visc0 * fac
       dh = dh0 * fac
       dv = dv0 * fac
       dc = dc0 * fac

       do ic = 1, ncan(p)

          if (dpai(p,ic) > 0._r8) then

             ! Reynolds number, Prandtl number, and Grashof number

             re = wind(p,ic) * dleaf(patch%itype(p)) / visc
             pr  = visc / dh
             gr = grav * dleaf(patch%itype(p))**3 * max(tleaf(p,ic,il)-tair(p,ic), 0._r8) / (tair(p,ic) * visc * visc)

             ! Nusselt number depends on convection regime

             ! a. Forced convection
             ! Note use of minimum conductance applied to gbh (needed for low wind speed)

             ! (i) Laminar flow

                nu = gb_factor * 0.66_r8 * pr**0.33_r8 * re**0.5_r8
                gbh_lam = (dh * nu / dleaf(patch%itype(p))) * rhomol(p) ; gbh_lam = max(gbh_lam, gbh_min)
                gbv_lam = gbh_lam * (dv / dh)**0.67_r8
                gbc_lam = gbh_lam * (dc / dh)**0.67_r8

             ! (ii) Turbulent flow

                nu = gb_factor * 0.036_r8 * pr**0.33_r8 * re**0.8_r8
                gbh_turb = (dh * nu / dleaf(patch%itype(p))) * rhomol(p) ; gbh_turb = max(gbh_turb, gbh_min)
                gbv_turb = gbh_turb * (dv / dh)**0.67_r8
                gbc_turb = gbh_turb * (dc / dh)**0.67_r8

             ! b. Free convection

                nu = 0.54_r8 * pr**0.25_r8 * gr**0.25_r8
                gbh_free = (dh * nu / dleaf(patch%itype(p))) * rhomol(p)
                gbv_free = gbh_free * (dv / dh)**0.75_r8
                gbc_free = gbh_free * (dc / dh)**0.75_r8

             ! Choose flow regimes to use

             select case (gb_type)
             case (1)

                ! Use only laminar flow

                gbh(p,ic,il) = gbh_lam
                gbv(p,ic,il) = gbv_lam
                gbc(p,ic,il) = gbc_lam

             case (2)

                ! Use laminar and turbulent flow

                gbh(p,ic,il) = max(gbh_lam, gbh_turb)
                gbv(p,ic,il) = max(gbv_lam, gbv_turb)
                gbc(p,ic,il) = max(gbc_lam, gbc_turb)

             case (3)

                ! Both forced and free convection occur together

                gbh(p,ic,il) = max(gbh_lam, gbh_turb) + gbh_free
                gbv(p,ic,il) = max(gbv_lam, gbv_turb) + gbv_free
                gbc(p,ic,il) = max(gbc_lam, gbc_turb) + gbc_free

             case default

                call endrun (msg=' ERROR: LeafBoundaryLayer: gb_type not valid')

             end select

          else

             gbh(p,ic,il) = 0._r8
             gbv(p,ic,il) = 0._r8
             gbc(p,ic,il) = 0._r8

          end if

       end do
    end do

    end associate
  end subroutine LeafBoundaryLayer

  !-----------------------------------------------------------------------
  subroutine CalcLeafBoundaryLayer (d, u, tleaf, tair, tref, pref, rhomol, gb_type_in, gbh, gbv, gbc)
    !
    ! !DESCRIPTION:
    ! Scalar leaf boundary layer conductance helper for testing.
    ! Mirrors the physics of LeafBoundaryLayer for a single leaf.
    !
    ! !USES:
    use clm_varcon, only : tfrz, grav
    use MLclm_varcon, only : visc0, dh0, dv0, dc0, gb_factor, gbh_min
    !
    ! !ARGUMENTS:
    implicit none
    real(r8), intent(in)  :: d         ! Leaf dimension (m)
    real(r8), intent(in)  :: u         ! Wind speed (m/s)
    real(r8), intent(in)  :: tleaf     ! Leaf temperature (K)
    real(r8), intent(in)  :: tair      ! Air temperature at leaf (K)
    real(r8), intent(in)  :: tref      ! Air temperature at reference height (K)
    real(r8), intent(in)  :: pref      ! Air pressure at reference height (Pa)
    real(r8), intent(in)  :: rhomol    ! Molar density of air (mol/m3)
    integer,  intent(in)  :: gb_type_in  ! Convection regime selector (1/2/3)
    real(r8), intent(out) :: gbh       ! Boundary layer conductance: heat (mol/m2/s)
    real(r8), intent(out) :: gbv       ! Boundary layer conductance: H2O (mol/m2/s)
    real(r8), intent(out) :: gbc       ! Boundary layer conductance: CO2 (mol/m2/s)
    !
    ! !LOCAL VARIABLES:
    real(r8) :: visc, dh, dv, dc, fac
    real(r8) :: nu, pr, re, gr
    real(r8) :: gbh_lam, gbv_lam, gbc_lam
    real(r8) :: gbh_turb, gbv_turb, gbc_turb
    real(r8) :: gbh_free, gbv_free, gbc_free
    !---------------------------------------------------------------------

    fac  = 101325._r8 / pref * (tref / tfrz)**1.81_r8
    visc = visc0 * fac
    dh   = dh0   * fac
    dv   = dv0   * fac
    dc   = dc0   * fac

    re = u * d / visc
    pr = visc / dh
    gr = grav * d**3 * max(tleaf - tair, 0._r8) / (tair * visc * visc)

    ! Laminar forced convection
    nu       = gb_factor * 0.66_r8 * pr**0.33_r8 * re**0.5_r8
    gbh_lam  = (dh * nu / d) * rhomol ; gbh_lam  = max(gbh_lam,  gbh_min)
    gbv_lam  = gbh_lam  * (dv / dh)**0.67_r8
    gbc_lam  = gbh_lam  * (dc / dh)**0.67_r8

    ! Turbulent forced convection
    nu       = gb_factor * 0.036_r8 * pr**0.33_r8 * re**0.8_r8
    gbh_turb = (dh * nu / d) * rhomol ; gbh_turb = max(gbh_turb, gbh_min)
    gbv_turb = gbh_turb * (dv / dh)**0.67_r8
    gbc_turb = gbh_turb * (dc / dh)**0.67_r8

    ! Free convection
    nu       = 0.54_r8 * pr**0.25_r8 * gr**0.25_r8
    gbh_free = (dh * nu / d) * rhomol
    gbv_free = gbh_free * (dv / dh)**0.75_r8
    gbc_free = gbh_free * (dc / dh)**0.75_r8

    select case (gb_type_in)
    case (1)
       gbh = gbh_lam  ; gbv = gbv_lam  ; gbc = gbc_lam
    case (2)
       gbh = max(gbh_lam, gbh_turb)
       gbv = max(gbv_lam, gbv_turb)
       gbc = max(gbc_lam, gbc_turb)
    case (3)
       gbh = max(gbh_lam, gbh_turb) + gbh_free
       gbv = max(gbv_lam, gbv_turb) + gbv_free
       gbc = max(gbc_lam, gbc_turb) + gbc_free
    case default
       call endrun (msg=' ERROR: CalcLeafBoundaryLayer: gb_type_in not valid')
    end select

  end subroutine CalcLeafBoundaryLayer

end module MLLeafBoundaryLayerMod
