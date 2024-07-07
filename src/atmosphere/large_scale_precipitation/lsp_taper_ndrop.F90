! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Large-scale precipitation scheme. Cloud droplet number calculator
! Subroutine Interface:
MODULE lsp_taper_ndrop_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='LSP_TAPER_NDROP_MOD'

CONTAINS


SUBROUTINE lsp_taper_ndrop(                                                   &
                   ! (Full) Aerosol tracers
                            so4_ait, so4_acc, so4_dis, sea_salt_film,         &
                            biogenic, sea_salt_jet, bmass_agd,                &
                            bmass_cld, ocff_agd, ocff_cld,                    &
                            nitr_acc, nitr_diss,                              &
                            n_arcl_species, n_arcl_compnts,                   &
                            i_arcl_compnts, arcl,                             &
                   ! Murk aerosol
                            aerosol,                                          &
                   ! CDNC from UKCA
                            ukca_cdnc,                                        &
                            cdnc_dim1, cdnc_dim2, cdnc_dim3,                  &
                   ! CDNC from EasyAerosol
                            easyaerosol_cdnc,                                 &
                   ! Other parameters
                            rhodz_dry, rhodz_moist,                           &
                            deltaz,                                           &
                            snow_depth, land_fract,                           &
                   ! Output parameter of n_drop_tpr
                            n_drop_tpr                                        &
                                 )
! Microphysics modules

USE mphys_inputs_mod,      ONLY: ndrop_surf,                                  &
                                 l_autoconv_murk, l_taper_new,                &
                                 max_drop_surf, l_droplet_tpr,                &
                                 l_mcr_arcl, arcl_inhom_sc,                   &
                                 l_use_sulphate_autoconv,                     &
                                 l_use_bmass_autoconv,                        &
                                 l_use_ocff_autoconv,                         &
                                 l_use_nitrate_autoconv,                      &
                                 l_use_seasalt_autoconv

USE mphys_constants_mod,   ONLY: ntot_land, ntot_sea, max_drop,               &
                                 n0_murk, m0_murk

USE lsp_autoc_consts_mod,  ONLY: power_murk, min_drop_alt, eta_peak,          &
                                 eta_low_nd, level_peak, level_surf,          &
                                 vala_fac1, vala_fac2,                        &
                                 half_range

! Dynamics and grid bounds modules

USE level_heights_mod,     ONLY: eta_theta_levels
USE atm_fields_bounds_mod, ONLY: tdims

! General constants modules
USE arcl_mod,              ONLY: npd_arcl_compnts, ip_arcl_sulp_ak,           &
                                 ip_arcl_sulp_ac,  ip_arcl_sulp_di,           &
                                 ip_arcl_sslt_fi,  ip_arcl_sslt_jt,           &
                                 ip_arcl_biom_ag, ip_arcl_biom_ic,            &
                                 ip_arcl_ocff_ag, ip_arcl_ocff_ic

USE conversions_mod,       ONLY: pi
USE rad_input_mod,         ONLY: l_use_biogenic, l_use_arclbiom,              &
                                 l_use_arclsslt, l_use_arclocff
USE gen_phys_inputs_mod,   ONLY: l_mr_physics
USE ukca_option_mod,       ONLY: l_ukca_aie2
USE glomap_clim_option_mod, ONLY: l_glomap_clim_aie2
USE murk_inputs_mod,       ONLY: l_murk
USE number_droplet_mod,    ONLY: number_droplet

! Stochastic physics module
USE stochastic_physics_run_mod, ONLY: l_rp2, i_rp_scheme, i_rp2b,             &
                                      ndrop_surf_rp

! EasyAerosol
USE easyaerosol_option_mod,ONLY: l_easyaerosol_autoconv
USE def_easyaerosol, ONLY: t_easyaerosol_cdnc

! Dr Hook modules
!--------------------------------------
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
!--------------------------------------

IMPLICIT NONE

!---------------------------------------------------------------------
! Purpose:
!   Calculates cloud drop number concentration for the microphysics
!   scheme.

!   By using the tapering method, we can also reduce droplet number
!   in the atmospheric boundary layer, from a peak at a given altitude
!   (z_peak_nd) to a user defined value (ndrop_surf) or a variable
!   value dependent on aerosol amounts.

! Method:
!  1) Determine the height (specifically eta value) below which droplet
!     number is to taper. When the taper curve is inactive, this value
!     is not used.

!  2) Above this height, calculate droplet number using the
!     Haywood-Jones formulae for Murk aerosol
!     (as in the autoconversion routine).

!     Haywood et al (2008): (aerosol number) :
!              n_aer = n0_murk * ( Aerosol / m0_murk ) ** power_murk

!     Jones  et al (1994): Aerosol number to droplet number:
!              n_d   = 3.75e8 * ( 1 - exp ( -2.5e-9 * n_aer ) )

!     If full prognostic aerosols are used or aerosol climatologies
!     are used, then the routine calls the number_droplet function
!     to generate cloud droplet number concentration.

!     If no aerosol species are available and droplet tapering is
!     requested, then the routine uses a simple profile of
!     pre-determined values for droplet number.

!  3) Using the first height above z_peak_nd, determine the droplet
!     profile at each level towards the surface

!    n_drop = ndrop_surf + vala * log ( z / z_surf )

!    where vala = ( nd(z_peak_nd) - ndrop_surf ) / log ( z_peak_nd / z_surf)

!    z_surf is the height at which the cloud drop number reaches ndrop_surf
!    and it is kept constant below this

!    ndrop_surf may be fixed (from user input) or be calculated as a
!    variable value if l_taper_new is set to .TRUE.. In this latter case,
!    n_drop = ndrop_surf2 + vala * log ( z / z_surf )

!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Large Scale Precipitation

! Description of Code:
!   FORTRAN 77  + common extensions also in Fortran90.
!   This code is written to UMDP3 version 6 programming standards.

!   Documentation: UMDP 26.

!------------------------------------------------
! Subroutine arguments
!------------------------------------------------

INTEGER :: cdnc_dim1, cdnc_dim2, cdnc_dim3

INTEGER, INTENT(IN) ::                                                        &
  n_arcl_species,                                                             &
                    ! Number of requested species within the climatology
  n_arcl_compnts,                                                             &
                    ! Corresponding number of requested components
  i_arcl_compnts(npd_arcl_compnts)
                    ! Array index of each aerosol clim component

REAL, INTENT(IN) :: ukca_cdnc(cdnc_dim1, cdnc_dim2, cdnc_dim3)
! CDNC from UKCA

TYPE (t_easyaerosol_cdnc), INTENT(IN) :: easyaerosol_cdnc
! CDNC from EasyAerosol

REAL, INTENT(IN) ::                                                           &
!-----------
! Aerosols
!-----------
  so4_ait      ( tdims%i_start : tdims%i_end,                                 &
                 tdims%j_start : tdims%j_end,                                 &
                             1 : tdims%k_end ),                               &
! SO4
! Aitken
  so4_acc      ( tdims%i_start : tdims%i_end,                                 &
                 tdims%j_start : tdims%j_end,                                 &
                             1 : tdims%k_end ),                               &
! SO4
! acc
  so4_dis      ( tdims%i_start : tdims%i_end,                                 &
                 tdims%j_start : tdims%j_end,                                 &
                             1 : tdims%k_end ),                               &
! SO4
! dis
  sea_salt_film( tdims%i_start : tdims%i_end,                                 &
                 tdims%j_start : tdims%j_end,                                 &
                             1 : tdims%k_end ),                               &
! Sea salt
! film
  biogenic     ( tdims%i_start : tdims%i_end,                                 &
                 tdims%j_start : tdims%j_end,                                 &
                             1 : tdims%k_end ),                               &
! Biogenic
! aerosol
  sea_salt_jet ( tdims%i_start : tdims%i_end,                                 &
                 tdims%j_start : tdims%j_end,                                 &
                             1 : tdims%k_end ),                               &
! Sea salt
! Jet
  bmass_agd    ( tdims%i_start : tdims%i_end,                                 &
                 tdims%j_start : tdims%j_end,                                 &
                             1 : tdims%k_end ),                               &
! Aged
! Biomass
  bmass_cld    ( tdims%i_start : tdims%i_end,                                 &
                 tdims%j_start : tdims%j_end,                                 &
                             1 : tdims%k_end ),                               &
! Cloudy
! Biomass
  ocff_agd     ( tdims%i_start : tdims%i_end,                                 &
                 tdims%j_start : tdims%j_end,                                 &
                             1 : tdims%k_end ),                               &
! ocff
! aged
  ocff_cld     ( tdims%i_start : tdims%i_end,                                 &
                 tdims%j_start : tdims%j_end,                                 &
                             1 : tdims%k_end ),                               &
! ocff
! cloud
  nitr_acc     ( tdims%i_start : tdims%i_end,                                 &
                 tdims%j_start : tdims%j_end,                                 &
                             1 : tdims%k_end ),                               &
! Nitr
! aged
  nitr_diss    ( tdims%i_start : tdims%i_end,                                 &
                 tdims%j_start : tdims%j_end,                                 &
                             1 : tdims%k_end ),                               &
!nitrate
!diss
  aerosol      ( tdims%i_start : tdims%i_end,                                 &
                 tdims%j_start : tdims%j_end,                                 &
                             1 : tdims%k_end ),                               &
! Murk aerosol for
! droplet number
! calculations
  arcl         ( tdims%i_start : tdims%i_end,                                 &
                 tdims%j_start : tdims%j_end,                                 &
                             1 : tdims%k_end,                                 &
                             1 : n_arcl_compnts ),                            &
! Aerosol climatologies for droplet number calculations

  rhodz_dry    ( tdims%i_start : tdims%i_end,                                 &
                 tdims%j_start : tdims%j_end,                                 &
                             1 : tdims%k_end ),                               &
  rhodz_moist  ( tdims%i_start : tdims%i_end,                                 &
                 tdims%j_start : tdims%j_end,                                 &
                             1 : tdims%k_end ),                               &
  deltaz       ( tdims%i_start : tdims%i_end,                                 &
                 tdims%j_start : tdims%j_end,                                 &
                             1 : tdims%k_end ),                               &
!Grid box information

  land_fract   ( tdims%i_start : tdims%i_end,                                 &
                 tdims%j_start : tdims%j_end ),                               &
! Land fraction
!(for Jones et al,
! 1994 method)
  snow_depth   ( tdims%i_start : tdims%i_end ,                                &
                 tdims%j_start : tdims%j_end )
! Depth of snow
!(for Jones et al,
! 1994 method)

REAL, INTENT(OUT) ::                                                          &
! intended output
    n_drop_tpr( tdims%i_start : tdims%i_end,                                  &
                tdims%j_start : tdims%j_end, 1:tdims%k_end )
! Droplet number after tapering has taken place

!------------------------------------------------
! Local variables
!------------------------------------------------

REAL ::                                                                       &
      n_aer( tdims%i_start : tdims%i_end,                                     &
             tdims%j_start : tdims%j_end  ),                                  &
         ! Aerosol number from
         ! Haywood et al (2008)
      vala,n_aer2,                                                            &
         ! Droplet number determined constants
      ndrop_surf2( tdims%i_start : tdims%i_end,                               &
                   tdims%j_start : tdims%j_end )
         ! Variable surface droplet number

! 3D version of Air density in kg/m3
REAL :: rho3d(  tdims%i_start : tdims%i_end,                                  &
                tdims%j_start : tdims%j_end,                                  &
                1 : tdims%k_end )

! A dummy array
REAL :: dummy( tdims%i_start : tdims%i_end,                                   &
               tdims%j_start : tdims%j_end,                                   &
               1 : tdims%k_end)

INTEGER :: i, j, k              ! Loop counters
INTEGER :: vec_length

! Logical to set nitrate climatology. Currently hardwired to .false.
! as a nitrate climatology is not yet available.
LOGICAL, PARAMETER :: l_use_arclnitr = .FALSE.

! Is the input "sulphate" aerosol in the form of
! ammonium sulphate (T) or just sulphur (F)? Used in the
! same way for the nitrate aerosol, in form of ammonium
! nitrate (T) or just nitrogen (F).
LOGICAL, PARAMETER :: l_nh42so4 = .FALSE.

REAL :: temp1, temp2
! temp for plane invariant reciprocals
REAL ::  tempv  
! temp  variable for intermediate computations

!------------------------------------------------
! Dr Hook subroutine timer details
!------------------------------------------------
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='LSP_TAPER_NDROP'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
!---------------------------------------------------------------------
! Start of physics
!---------------------------------------------------------------------

! First update value of ndrop_surf if random parameters 2b
! scheme is switched on
IF ( l_rp2 .AND. i_rp_scheme == i_rp2b ) THEN
  ndrop_surf = ndrop_surf_rp
END IF

vec_length  = tdims%i_len*tdims%j_len

IF (l_droplet_tpr) THEN

  !------------------------------------------------
  ! Droplet tapering is on, so we need to calculate
  ! droplet number concentration and add in a taper
  !------------------------------------------------

  IF (l_murk .AND. l_autoconv_murk) THEN

    !----------------------------------------------
    ! Use murk aerosol to calculate droplet number
    ! and taper this profile
    !----------------------------------------------

    ! Calculate variable surface droplet number
    ! if required

    IF (l_taper_new) THEN

      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end

            ! Use Haywood et al (2008) formulae to get droplet number

          n_aer(i,j) = MAX( aerosol(i,j,1) / m0_murk*1.0e-9, 0.0001)
          ! 1.0E-9 converts from ug/kg to kg/kg

          !-----------------------------------------------
          ! Calculation of the aerosol number
          !-----------------------------------------------
          n_aer(i,j) = n0_murk * n_aer(i,j) ** power_murk

          !-----------------------------------------------
          ! Convert to CCN using a Jones et al (1994)
          ! relationship (as modified by Jonathan Wilkinson)
          !-----------------------------------------------

          ndrop_surf2(i,j) =                                                  &
             max_drop_surf * (1.0 - EXP( - 1.5e-9 * n_aer(i,j) ) )

          !------------------------------------------------
          ! Ensure the surface droplet number doesn't get
          ! below the minimum value
          !------------------------------------------------

          ndrop_surf2(i,j) = MAX( ndrop_surf2(i,j), ndrop_surf)


        END DO ! tdims%i
      END DO   ! tdims%j

    ELSE

!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(i,j) &
!$OMP SHARED(tdims,ndrop_surf2,ndrop_surf)
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          ndrop_surf2(i,j) = ndrop_surf
        END DO ! tdims%i
      END DO   ! tdims%j
!$OMP END PARALLEL DO

    END IF ! l_taper_new (new version of taper code)


!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(i,j,k,n_aer2,tempv) &
!$OMP SHARED(tdims,level_peak,aerosol,land_fract,snow_depth,n_drop_tpr, &
!$OMP  n0_murk,m0_murk)
    DO k = level_peak,tdims%k_end

      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end

          ! Above the taper level, so use Haywood-Jones formulae:

          n_aer2 = MAX( aerosol(i,j, k) / m0_murk * 1.0e-9, 0.0001)
          ! 1.0E-9 converts from ug/kg to kg/kg

      !-----------------------------------------------
      ! Calculation of the aerosol number
      !-----------------------------------------------

          tempv=n_aer2**power_murk

      !-----------------------------------------------
      ! Convert to CCN using the Jones et al (1994)
      ! relationship (follows number_droplet routine
      ! but reproduced here to avoid compiler issues).
      !-----------------------------------------------

          n_aer2 = -2.5e-9 *  n0_murk  * tempv

          n_drop_tpr(i,j,k) = max_drop * ( 1.0e+00 - EXP(n_aer2)  )

          IF ( land_fract(i,j) >= 0.2 .AND.                                   &
               snow_depth(i,j) < 5000.0    ) THEN

            n_drop_tpr(i,j,k) = MAX( n_drop_tpr(i,j,k), 35.0e+06 )

          ELSE

            n_drop_tpr(i,j,k) = MAX( n_drop_tpr(i,j,k), 5.0e+06  )

          END IF ! land fract > 0.2 etc

        END DO
      END DO
    END DO
!$OMP END PARALLEL DO

    ! Start tapering

    DO k=1, level_peak-1

      temp1 =  LOG( eta_theta_levels(k) / eta_theta_levels(1) )
      temp2 =  LOG( eta_theta_levels(k) / eta_theta_levels(level_surf) )

      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end

          IF (n_drop_tpr(i,j,level_peak) > ndrop_surf2(i,j)) THEN

            ! if drop number is increasing with height, keep it at
            ! its surface value below z_surf before tapering up to
            ! the z_peak value

            vala = (n_drop_tpr(i,j,level_peak) - ndrop_surf2(i,j))            &
                   * vala_fac2

            n_drop_tpr( i, j, k ) =  MAX( ndrop_surf2(i,j),                   &
                 ndrop_surf2(i,j) + vala * temp2 )

          ELSE

            ! if drop number is decreasing with height, make sure it
            ! only reaches its surface value at level 1

            vala = (n_drop_tpr(i,j,level_peak) - ndrop_surf2(i,j))            &
                   * vala_fac1

            n_drop_tpr( i, j, k ) =  ndrop_surf2(i,j) + vala * temp1

          END IF

        END DO ! tdims%i
      END DO   ! tdims%j

    END DO ! eta values below peak

  ELSE IF (l_mcr_arcl) THEN

    !------------------------------------------------------------
    ! Full droplet number calculation with aerosol climatologies
    !------------------------------------------------------------
    ! Step 1: Calculate air density, rho
!$OMP  PARALLEL DEFAULT(NONE)                                          &
!$OMP& SHARED( level_peak, tdims, l_mr_physics, rho3d, rhodz_dry,     &
!$OMP&         deltaz, rhodz_moist )                                   &
!$OMP& PRIVATE( i, j, k )
    IF (l_mr_physics) THEN
!$OMP DO SCHEDULE(STATIC)
      DO k = level_peak, tdims%k_end
        DO j = tdims%j_start, tdims%j_end
          DO i = tdims%i_start, tdims%i_end
            ! rho is the dry density
            rho3d(i,j,k) = rhodz_dry(i,j,k) / deltaz(i,j,k)
          END DO
        END DO
      END DO
!$OMP END DO

    ELSE ! l_mr_physics
!$OMP DO SCHEDULE(STATIC)
      DO k = level_peak, tdims%k_end
        DO j = tdims%j_start, tdims%j_end
          DO i = tdims%i_start, tdims%i_end
            ! rho is the moist density
            rho3d(i,j,k) = rhodz_moist(i,j,k) / deltaz(i,j,k)
          END DO
        END DO
      END DO
!$OMP END DO
    END IF  ! l_mr_physics
!$OMP END PARALLEL

    ! Step 2: Call number_droplet routine
    CALL number_droplet(                                                      &
                         tdims%i_start, tdims%i_end,                          &
                         tdims%j_start, tdims%j_end,                          &
                         1, tdims%k_end,                                      &
                         level_peak, tdims%k_end,                             &
                         l_mcr_arcl,                                          &
                         l_nh42so4,                                           &
                         arcl(:,:,:,i_arcl_compnts(ip_arcl_sulp_ac)),         &
                         arcl(:,:,:,i_arcl_compnts(ip_arcl_sulp_di)),         &
                         l_use_arclsslt,                                      &
                         arcl(:,:,:,i_arcl_compnts(ip_arcl_sslt_fi)),         &
                         arcl(:,:,:,i_arcl_compnts(ip_arcl_sslt_jt)),         &
                         l_use_biogenic, biogenic,                            &
                         l_use_arclbiom,                                      &
                         arcl(:,:,:,i_arcl_compnts(ip_arcl_biom_ag)),         &
                         arcl(:,:,:,i_arcl_compnts(ip_arcl_biom_ic)),         &
                         l_use_arclocff,                                      &
                         arcl(:,:,:,i_arcl_compnts(ip_arcl_ocff_ag)),         &
                         arcl(:,:,:,i_arcl_compnts(ip_arcl_ocff_ic)),         &
                         l_use_arclnitr,                                      &
                         dummy,                                               &
                         dummy,                                               &
                         rho3d,                                               &
                         snow_depth, land_fract,                              &
                         n_drop_tpr)

!$OMP  PARALLEL DEFAULT(NONE)                                          &
!$OMP& SHARED( level_peak, tdims, n_drop_tpr, arcl_inhom_sc,           &
!$OMP&         ndrop_surf, eta_theta_levels, level_surf, vala_fac1,    &
!$OMP&         vala_fac2 )                                             &
!$OMP& PRIVATE( i, j, k, temp1, temp2, vala )
!$OMP  DO SCHEDULE(STATIC)
    DO k = level_peak, tdims%k_end
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          ! scaling the n_drop by the inhomogeneity scaling:
          n_drop_tpr(i,j,k) = n_drop_tpr(i,j,k) * arcl_inhom_sc

        END DO
      END DO
    END DO
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE(STATIC)
    DO k = 1, level_peak-1
      temp1 =  LOG( eta_theta_levels(k) / eta_theta_levels(1) )
      temp2 =  LOG( eta_theta_levels(k) / eta_theta_levels(level_surf) )
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          IF (n_drop_tpr(i,j,level_peak) > ndrop_surf) THEN

            ! if drop number is increasing with height, keep it at
            ! its surface value below z_surf before tapering up to
            ! the z_peak value

            vala = (n_drop_tpr(i,j,level_peak) - ndrop_surf)                  &
                   * vala_fac2

            n_drop_tpr( i, j, k ) =  MAX( ndrop_surf,                         &
                 ndrop_surf + vala * temp2 )

          ELSE

            ! if drop number is decreasing with height, make sure it
            ! only reaches its surface value at level 1

            vala = (n_drop_tpr(i,j,level_peak) - ndrop_surf)                  &
                   * vala_fac1

            n_drop_tpr( i, j, k ) =  ndrop_surf + vala * temp1

          END IF
        END DO ! tdims%i
      END DO ! tdims%j
    END DO ! tdims%k
!$OMP END DO
!$OMP END PARALLEL

  ELSE IF (l_use_sulphate_autoconv .OR. l_ukca_aie2 .OR. &
           l_glomap_clim_aie2 .OR. l_easyaerosol_autoconv) THEN

    !------------------------------------------------------------
    ! Full droplet number calculation based on prognostic aerosol
    !------------------------------------------------------------

!$OMP  PARALLEL DEFAULT(NONE)                                          &
!$OMP& SHARED( l_ukca_aie2, l_glomap_clim_aie2, l_mr_physics, tdims,   &
!$OMP&         n_drop_tpr, rho3d,                                      &
!$OMP&         level_peak, ukca_cdnc, rhodz_dry, deltaz, rhodz_moist,  &
!$OMP&         l_easyaerosol_autoconv, easyaerosol_cdnc )              &
!$OMP& PRIVATE( i, j, k )
    IF (l_easyaerosol_autoconv) THEN
!$OMP DO SCHEDULE(STATIC)
      DO k = level_peak, tdims%k_end
        DO j = tdims%j_start, tdims%j_end
          DO i = tdims%i_start, tdims%i_end
            n_drop_tpr(i,j,k) = easyaerosol_cdnc%cdnc(i,j,k)
          END DO
        END DO
      END DO
!$OMP END DO
    ELSE
      IF (l_ukca_aie2 .OR. l_glomap_clim_aie2) THEN
!$OMP DO SCHEDULE(STATIC)
        DO k = level_peak, tdims%k_end
          DO j = tdims%j_start, tdims%j_end
            DO i = tdims%i_start, tdims%i_end
              n_drop_tpr(i,j,k) = ukca_cdnc(i,j,k)
            END DO
          END DO
        END DO
!$OMP END DO
      ELSE  ! (not l_ukca_aie2 AND not l_glomap_clim_aie2)

            ! Step 1: Calculate air density, rho

        IF (l_mr_physics) THEN
!$OMP DO SCHEDULE(STATIC)
          DO k = level_peak, tdims%k_end
            DO j = tdims%j_start, tdims%j_end
              DO i = tdims%i_start, tdims%i_end

                ! rho is the dry density

                rho3d(i,j,k) = rhodz_dry(i,j,k) / deltaz(i,j,k)
              END DO
            END DO
          END DO
!$OMP END DO
        ELSE ! l_mr_physics
!$OMP DO SCHEDULE(STATIC)
          DO k = level_peak, tdims%k_end
            DO j = tdims%j_start, tdims%j_end
              DO i = tdims%i_start, tdims%i_end

                ! rho is the moist density

                rho3d(i,j,k) = rhodz_moist(i,j,k) / deltaz(i,j,k)
              END DO
            END DO
          END DO
!$OMP END DO
        END IF  ! l_mr_physics

      END IF  ! l_ukca_aie2 .OR. l_glomap_clim_aie2

    END IF ! l_easyaerosol_autoconv
!$OMP END PARALLEL

    ! Step 2: Call number_droplet routine
    IF (.NOT. l_easyaerosol_autoconv .AND. .NOT. &
                                    (l_ukca_aie2 .OR. l_glomap_clim_aie2)) THEN

      CALL number_droplet(                                                    &
                       tdims%i_start, tdims%i_end,                            &
                       tdims%j_start, tdims%j_end,                            &
                       1, tdims%k_end,                                        &
                       level_peak, tdims%k_end,                               &
                       l_use_sulphate_autoconv,                               &
                       l_nh42so4,                                             &
                       so4_acc,                                               &
                       so4_dis,                                               &
                       l_use_seasalt_autoconv,                                &
                       sea_salt_film,                                         &
                       sea_salt_jet,                                          &
                       l_use_biogenic, biogenic,                              &
                       l_use_bmass_autoconv,                                  &
                       bmass_agd,                                             &
                       bmass_cld,                                             &
                       l_use_ocff_autoconv,                                   &
                       ocff_agd, ocff_cld,                                    &
                       l_use_nitrate_autoconv,                                &
                       nitr_acc,                                              &
                       nitr_diss,                                             &
                       rho3d,                                                 &
                       snow_depth, land_fract,                                &
                       n_drop_tpr)
    END IF


!$OMP  PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE)                      &
!$OMP& SHARED( level_peak, tdims, n_drop_tpr, ndrop_surf, level_surf,  &
!$OMP&         eta_theta_levels, vala_fac1, vala_fac2 )                &
!$OMP& PRIVATE( i, j, k, vala, temp1, temp2 )
    DO k = 1, level_peak-1
      temp1 =  LOG( eta_theta_levels(k) / eta_theta_levels(1) )
      temp2 =  LOG( eta_theta_levels(k) / eta_theta_levels(level_surf) )
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          IF (n_drop_tpr(i,j,level_peak) > ndrop_surf) THEN

            ! if drop number is increasing with height, keep it at
            ! its surface value below z_surf before tapering up to
            ! the z_peak value

            vala = (n_drop_tpr(i,j,level_peak) - ndrop_surf)                  &
                   * vala_fac2

            n_drop_tpr( i, j, k ) =  MAX( ndrop_surf,                         &
                 ndrop_surf + vala * temp2 )

          ELSE

            ! if drop number is decreasing with height, make sure it
            ! only reaches its surface value at level 1

            vala = (n_drop_tpr(i,j,level_peak) - ndrop_surf)                  &
                   * vala_fac1

            n_drop_tpr( i, j, k ) =  ndrop_surf + vala * temp1

          END IF
        END DO ! tdims%i
      END DO ! tdims%j
    END DO ! tdims%k
!$OMP END PARALLEL DO

  ELSE   ! not l_murk, l_mcr_arcl, l_use_sulphate_autoconv or l_ukca_aie2
         ! or l_glomap_clim_aie2 or l_easyaerosol_autoconv

    !-----------------------------------------
    ! Use a simple profile for tapering
    !-----------------------------------------

    ! Do not allow a peak value of eta to exceed that of the
    ! pre-determined low droplet number (usually around 2km)

    IF (eta_peak > eta_low_nd) eta_peak = eta_low_nd

!$OMP  PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE)                    &
!$OMP& SHARED( tdims, eta_theta_levels, eta_low_nd, eta_peak,        &
!$OMP&         n_drop_tpr, ndrop_surf, vala_fac1, vala_fac2,         &
!$OMP&         level_surf )                                          &
!$OMP& PRIVATE( i, j, k, vala, temp1, temp2 )
    DO k = 1, tdims%k_end

      IF ( eta_theta_levels(k) >= eta_low_nd .AND.                            &
           eta_theta_levels(k) >= eta_peak         ) THEN

        ! Above taper level and level set for minimum droplet
        ! number, hence set n_drop at this level to low value

        DO j = tdims%j_start, tdims%j_end
          DO i = tdims%i_start, tdims%i_end
            n_drop_tpr(i,j,k) = min_drop_alt
          END DO
        END DO

      ELSE IF ( eta_theta_levels(k) >= eta_peak .AND.                         &
                eta_theta_levels(k) < eta_low_nd     ) THEN

        ! Above taper level yet below level set for minimum
        ! droplet number, so use a simple function of eta to
        ! determine the droplet number. This is intended to
        ! be independent of model levels

        DO j = tdims%j_start, tdims%j_end
          DO i = tdims%i_start, tdims%i_end
            n_drop_tpr(i,j,k) = half_range * (- COS (pi + pi *                &
                          ( ( eta_theta_levels (k) - eta_peak ) /             &
                          ( eta_low_nd - eta_peak ) ) ) )                     &
                          + half_range + min_drop_alt
          END DO
        END DO

        ! No need to set up peak droplet as this should be max_drop

      ELSE ! eta_theta_levels

        IF (max_drop > ndrop_surf) THEN

          ! if drop number is increasing with height, keep it at
          ! its surface value below z_surf before tapering up to
          ! the z_peak value
          temp2 =  LOG( eta_theta_levels(k) / eta_theta_levels(level_surf) )
          vala = (max_drop - ndrop_surf) / vala_fac2

          DO j = tdims%j_start, tdims%j_end
            DO i = tdims%i_start, tdims%i_end
              n_drop_tpr( i, j, k ) =  MAX( ndrop_surf,                       &
                   ndrop_surf + vala * temp2 )
            END DO ! tdims%i
          END DO ! tdims%j

        ELSE

          ! if drop number is decreasing with height, make sure it
          ! only reaches its surface value at level 1
          temp1 =  LOG( eta_theta_levels(k) / eta_theta_levels(1) )
          vala = (max_drop - ndrop_surf) * vala_fac1

          DO j = tdims%j_start, tdims%j_end
            DO i = tdims%i_start, tdims%i_end
              n_drop_tpr( i, j, k ) =  ndrop_surf + vala * temp1
            END DO ! tdims%i
          END DO ! tdims%j

        END IF

      END IF !eta_levels above or below peak values

    END DO  ! tdims%k
!$OMP END PARALLEL DO

  END IF ! l_murk / l_use_sulphate_autoconv

ELSE ! l_droplet_tpr

  !-------------------------------------------------
  ! Droplet tapering is not active, but we still
  ! need to calculate potential cloud drop number
  ! concentration for use in autoconversion
  !-------------------------------------------------

  IF (l_murk .AND. l_autoconv_murk) THEN

    ! Calculate using MURK aerosol
    DO k = 1,tdims%k_end
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end

          ! Use Jones-Haywood Formulae
          n_aer(i,j) =                                                        &
                     MAX( aerosol(i,j, k) / m0_murk * 1.0e-9, 0.0001)
          ! 1.0E-9 converts from ug/kg to kg/kg

          !-----------------------------------------------
          ! Calculation of the aerosol number
          !-----------------------------------------------
          n_aer(i,j) = n0_murk * n_aer(i,j) ** power_murk
          n_drop_tpr(i,j,k) =                                                 &
              max_drop * ( 1.0e+00-EXP( -2.5e-9 * n_aer(i,j) ) )

          IF ( land_fract(i,j) >= 0.2 .AND.                                   &
               snow_depth(i,j) < 5000.0    ) THEN

            n_drop_tpr(i,j,k) = MAX( n_drop_tpr(i,j,k), 35.0e+06 )

          ELSE

            n_drop_tpr(i,j,k) = MAX( n_drop_tpr(i,j,k), 5.0e+06  )

          END IF ! land fract > 0.2 etc

        END DO ! i (tdims%i)
      END DO ! j (tdims%j)
    END DO ! k (tdims%k)

  ELSE IF (l_mcr_arcl) THEN

    ! Calculate using aerosol climatologies

!$OMP  PARALLEL DEFAULT(NONE)                                        &
!$OMP& SHARED( l_mr_physics, tdims, rho3d, rhodz_dry, rhodz_moist,  &
!$OMP&         deltaz )                                              &
!$OMP& PRIVATE( i, j, k )
    IF (l_mr_physics) THEN
!$OMP DO SCHEDULE(STATIC)
      DO k = 1,tdims%k_end
        DO j = tdims%j_start, tdims%j_end
          DO i = tdims%i_start, tdims%i_end

          ! Step 1: Calculate air density, rho

            ! rho is the dry density
            rho3d(i,j,k) = rhodz_dry(i,j,k) / deltaz(i,j,k)
          END DO
        END DO
      END DO
!$OMP END DO

    ELSE ! l_mr_physics

!$OMP DO SCHEDULE(STATIC)
      DO k = 1,tdims%k_end
        DO j = tdims%j_start, tdims%j_end
          DO i = tdims%i_start, tdims%i_end

            ! rho is the moist density
            rho3d(i,j,k) = rhodz_moist(i,j,k) / deltaz(i,j,k)

          END DO
        END DO
      END DO
!$OMP END DO
    END IF  ! l_mr_physics
!$OMP END PARALLEL

    ! Step 2: Call number_droplet routine
    CALL number_droplet(                                                      &
                         tdims%i_start, tdims%i_end,                          &
                         tdims%j_start, tdims%j_end,                          &
                         1, tdims%k_end,                                      &
                         1, tdims%k_end,                                      &
                         l_mcr_arcl,                                          &
                         l_nh42so4,                                           &
                         arcl(:,:,:,i_arcl_compnts(ip_arcl_sulp_ac)),         &
                         arcl(:,:,:,i_arcl_compnts(ip_arcl_sulp_di)),         &
                         l_use_arclsslt,                                      &
                         arcl(:,:,:,i_arcl_compnts(ip_arcl_sslt_fi)),         &
                         arcl(:,:,:,i_arcl_compnts(ip_arcl_sslt_jt)),         &
                         l_use_biogenic, biogenic,                            &
                         l_use_arclbiom,                                      &
                         arcl(:,:,:,i_arcl_compnts(ip_arcl_biom_ag)),         &
                         arcl(:,:,:,i_arcl_compnts(ip_arcl_biom_ic)),         &
                         l_use_arclocff,                                      &
                         arcl(:,:,:,i_arcl_compnts(ip_arcl_ocff_ag)),         &
                         arcl(:,:,:,i_arcl_compnts(ip_arcl_ocff_ic)),         &
                         l_use_arclnitr,                                      &
                         dummy,                                               &
                         dummy,                                               &
                         rho3d,                                               &
                         snow_depth, land_fract,                              &
                         n_drop_tpr)

!$OMP  PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE)                    &
!$OMP& SHARED( tdims, n_drop_tpr, arcl_inhom_sc )                    &
!$OMP& PRIVATE( i, j, k )
    DO k = 1,tdims%k_end
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          ! scaling the n_drop by the inhomogeneity scaling:
          n_drop_tpr(i,j,k) = n_drop_tpr(i,j,k) * arcl_inhom_sc

        END DO
      END DO
    END DO
!$OMP END PARALLEL DO

  ELSE IF (l_use_sulphate_autoconv .OR. l_ukca_aie2 .OR. &
           l_glomap_clim_aie2 .OR. l_easyaerosol_autoconv) THEN

    !------------------------------------------------------------
    ! Full droplet number calculation based on prognostic aerosol
    !------------------------------------------------------------

!$OMP  PARALLEL DEFAULT(NONE)                                        &
!$OMP& SHARED( l_ukca_aie2, l_glomap_clim_aie2, l_mr_physics, tdims, &
!$OMP&         n_drop_tpr,                                           &
!$OMP&         ukca_cdnc, rho3d, rhodz_dry, deltaz, rhodz_moist,     &
!$OMP&         l_easyaerosol_autoconv, easyaerosol_cdnc )            &
!$OMP& PRIVATE( i, j, k )
    IF (l_easyaerosol_autoconv) THEN
!$OMP DO SCHEDULE(STATIC)
      DO k = 1, tdims%k_end
        DO j = tdims%j_start, tdims%j_end
          DO i = tdims%i_start, tdims%i_end
            n_drop_tpr(i,j,k) = easyaerosol_cdnc%cdnc(i,j,k)
          END DO
        END DO
      END DO
!$OMP END DO
    ELSE
      IF (l_ukca_aie2 .OR. l_glomap_clim_aie2) THEN
!$OMP DO SCHEDULE(STATIC)
        DO k = 1, tdims%k_end
          DO j = tdims%j_start, tdims%j_end
            DO i = tdims%i_start, tdims%i_end
              n_drop_tpr(i,j,k) = ukca_cdnc(i,j,k)
            END DO
          END DO
        END DO
!$OMP END DO
      ELSE  ! (not l_easyaerosol_autoconv and not 
            !                               (l_ukca_aie2 or l_glomap_clim_aie2)

            ! Step 1: Calculate air density, rho

        IF (l_mr_physics) THEN
!$OMP DO SCHEDULE(STATIC)
          DO k = 1, tdims%k_end
            DO j = tdims%j_start, tdims%j_end
              DO i = tdims%i_start, tdims%i_end

                ! rho is the dry density

                rho3d(i,j,k) = rhodz_dry(i,j,k) / deltaz(i,j,k)
              END DO
            END DO
          END DO
!$OMP END DO
        ELSE ! l_mr_physics
!$OMP DO SCHEDULE(STATIC)
          DO k = 1, tdims%k_end
            DO j = tdims%j_start, tdims%j_end
              DO i = tdims%i_start, tdims%i_end

                ! rho is the moist density

                rho3d(i,j,k) = rhodz_moist(i,j,k) / deltaz(i,j,k)
              END DO
            END DO
          END DO
!$OMP END DO
        END IF  ! l_mr_physics

      END IF  ! l_ukca_aie2 or l_glomap_clim_aie2

    END IF ! l_easyaerosol_autoconv
!$OMP END PARALLEL

    IF (.NOT. (l_ukca_aie2 .OR. l_glomap_clim_aie2) .AND. .NOT. &
                                                   l_easyaerosol_autoconv) THEN

      ! Step 2: Call number_droplet routine
      CALL number_droplet(                                                    &
                         tdims % i_start, tdims % i_end,                      &
                         tdims % j_start, tdims % j_end,                      &
                         1, tdims % k_end,                                    &
                         1, tdims % k_end,                                    &
                         l_use_sulphate_autoconv,                             &
                         l_nh42so4,                                           &
                         so4_acc,                                             &
                         so4_dis,                                             &
                         l_use_seasalt_autoconv,                              &
                         sea_salt_film,                                       &
                         sea_salt_jet,                                        &
                         l_use_biogenic, biogenic,                            &
                         l_use_bmass_autoconv,                                &
                         bmass_agd,                                           &
                         bmass_cld,                                           &
                         l_use_ocff_autoconv,                                 &
                         ocff_agd, ocff_cld,                                  &
                         l_use_nitrate_autoconv,                              &
                         nitr_acc,                                            &
                         nitr_diss,                                           &
                         rho3d,                                               &
                         snow_depth, land_fract,                              &
                         n_drop_tpr)
    END IF ! l_ukca_aie2 or l_glomap_clim_aie2

  ELSE ! (not l_murk, l_mcr_arcl, l_use_sulphate_autoconv or l_ukca_aie2 or
       !  l_glomap_clim_aie2 or l_easyaerosol_autoconv)

    !-----------------------------------------------------------
    ! In this case, STASH 4/211 will be a simple land-sea split
    ! dependant on ntot_land and ntot_sea
    !-----------------------------------------------------------

    ! Set the first level
!$OMP  PARALLEL DEFAULT(NONE)                                        &
!$OMP& SHARED( tdims, land_fract, n_drop_tpr )                       &
!$OMP& PRIVATE( i, j, k )
!$OMP DO SCHEDULE(STATIC)
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end

        IF (land_fract(i,j) >= 0.5) THEN
          n_drop_tpr(i,j,1) = ntot_land
        ELSE ! land_fract
          n_drop_tpr(i,j,1) = ntot_sea
        END IF ! land fract

      END DO ! i (tdims%i)
    END DO ! j (tdims%j)
!$OMP END DO

    ! copy up to higher levels - thus avoiding too many branches
    ! in loops
!$OMP DO SCHEDULE(STATIC)
    DO k = 2, tdims%k_end
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          n_drop_tpr(i,j,k) = n_drop_tpr(i,j,1)
        END DO
      END DO
    END DO
!$OMP END DO
!$OMP END PARALLEL

  END IF ! l_autoconv_murk, l_use_sulphate_autoconv

END IF ! l_droplet_tpr

!---------------------------------------------------------------------
! End of physics
!---------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE lsp_taper_ndrop
END MODULE lsp_taper_ndrop_mod
