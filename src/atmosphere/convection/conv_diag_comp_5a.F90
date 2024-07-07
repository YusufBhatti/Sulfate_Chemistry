! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! To diagnose convective occurrence and type
!
MODULE conv_diag_comp_5a_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'CONV_DIAG_COMP_5A_MOD'
CONTAINS

SUBROUTINE conv_diag_comp_5a(                                           &

! IN values defining field dimensions and subset to be processed :
       row_length, rows                                                 &

! IN level and points info
      , bl_levels,  nunstable, index_i, index_j                         &

! IN grid information
      , p, P_theta_lev, exner_rho                                       &
      , rho_only, rho_theta, z_full, z_half                             &

! IN Cloud data :
      , qcf, qcl, cloud_fraction                                        &

! IN everything not covered so far :

      , pstar, q, theta, exner_theta_levels                             &
      , land_mask, frac_land, timestep                                  &
      , w_copy, w_max, tv1_sd, bl_vscale2                               &
      , deep_flag, past_precip, past_conv_ht                            &

! SCM Diagnostics (dummy values in full UM)
      , nSCMDpkgs, L_SCMDiags                                           &

! INOUT data required elsewhere in UM system :

      , ntml, ntpar, nlcl                                               &
      , cumulus, L_shallow, l_congestus, l_congestus2                   &
      , zh,zhpar,dzh,qcl_inv_top,z_lcl,z_lcl_uv,delthvu,ql_ad           &
      , cin, cape, entrain_coef                                         &
      , qsat_lcl                                                        &
       )

! Modules used
! Definitions of prognostic variable array sizes
USE atm_fields_bounds_mod, ONLY:                                        &
  pdims_s, tdims_s, tdims, ScmRowLen, ScmRow

! Model level heights from centre of Earth
USE level_heights_mod, ONLY: &
  r_theta_levels             &  ! Radii on theta levels (m)
 ,r_rho_levels                  ! Radii on rho levels (m)

USE cv_derived_constants_mod, ONLY:                                     &
  ls, lsrcp, lcrcp, gamma_dry

USE cv_diag_param_mod, ONLY:                                            &
  a_plume, b_plume, max_t_grad

USE cv_param_mod, ONLY:                                                 &
  max_diag_thpert, sh_grey_closure

USE bl_option_mod, ONLY: bl_res_inv, on

USE planet_constants_mod, ONLY:                                         &
    cp, r, repsilon, c_virtual, g, planet_radius
USE water_constants_mod, ONLY: lc, lf, tm

USE cv_run_mod, ONLY:                                                   &
    iconv_congestus, iconv_deep, icvdiag, cldbase_opt_dp,               &
    w_cape_limit, limit_pert_opt, cvdiag_inv, cvdiag_sh_wtest,          &
    cldbase_opt_sh

USE cloud_inputs_mod, ONLY: forced_cu

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

! subroutines
USE cumulus_test_mod, ONLY: cumulus_test
USE mean_w_layer_mod, ONLY: mean_w_layer
USE parcel_ascent_mod, ONLY: parcel_ascent

USE UM_ParParams
USE missing_data_mod, ONLY: rmdi, imdi
USE ereport_mod,      ONLY: ereport
USE model_domain_mod, ONLY: model_type, mt_single_column
USE s_scmop_mod,      ONLY: default_streams,                            &
                            t_inst, d_wet, d_point,                     &
                            scmdiag_bl, scmdiag_conv
USE scmoutput_mod,    ONLY: scmoutput

USE nlsizes_namelist_mod,   ONLY: model_levels

USE errormessagelength_mod, ONLY: errormessagelength

USE cv_parcel_neutral_dil_mod, ONLY: cv_parcel_neutral_dil
USE cv_parcel_neutral_inv_mod, ONLY: cv_parcel_neutral_inv
USE lift_cond_lev_mod, ONLY: lift_cond_lev
IMPLICIT NONE

! ------------------------------------------------------------------------------
! Description:
!   To diagnose convection occurrence and type on just unstable points.
!
!   Called by conv_diag
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Convection
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.1.
!
! ------------------------------------------------------------------------------
! Subroutine arguments

! (a) Defining horizontal grid and subset thereof to be processed.

INTEGER, INTENT(IN) :: &
  row_length           & ! Local number of points on a row
 ,rows                   ! Local number of rows in a theta field

! (b) Defining vertical grid of model atmosphere.

INTEGER, INTENT(IN) :: &
  bl_levels            & ! Max. no. of "boundary" levels allowed.
 ,nunstable              ! number of unstable points

INTEGER, INTENT(IN) :: &
  index_i(nunstable)   & ! column number of unstable points
 ,index_j(nunstable)     ! row number of unstable points

REAL, INTENT(IN) ::                   &
  p(pdims_s%i_start:pdims_s%i_end,    & ! pressure  on rho levels (Pa)
    pdims_s%j_start:pdims_s%j_end,    &
    pdims_s%k_start:pdims_s%k_end)    &
 ,p_theta_lev(tdims%i_start:tdims%i_end,    & ! P on theta lev (Pa)
              tdims%j_start:tdims%j_end,    &
                          1:tdims%k_end)    &
 ,exner_rho(pdims_s%i_start:pdims_s%i_end,  & ! Exner on rho level
            pdims_s%j_start:pdims_s%j_end,  & !
            pdims_s%k_start:pdims_s%k_end)  &
 ,exner_theta_levels(tdims%i_start:tdims%i_end, & ! exner pressure theta lev
                     tdims%j_start:tdims%j_end, & !  (Pa)
                                 1:tdims%k_end)

REAL, INTENT(IN) ::                         &
  rho_only(row_length,rows,model_levels)    & ! density (kg/m3)
 ,rho_theta(row_length,rows,model_levels-1) & ! density th lev (kg/m3)
 ,z_full(row_length,rows,model_levels)      & ! height th lev (m)
 ,z_half(row_length,rows,model_levels)        ! height rho lev (m)

! (c) Cloud data.
REAL, INTENT(IN) ::                         &
  qcf(row_length,rows,model_levels)         & ! Cloud ice (kg/kg air)
 ,qcl(row_length,rows,model_levels)         & ! Cloud liquid water (kg/kg air)
 ,cloud_fraction(row_length, rows, model_levels)   ! Cloud fraction

! (d) Atmospheric + any other data not covered so far, incl control.
REAL, INTENT(IN) ::                      &
  pstar(row_length, rows)                & ! Surface pressure (Pascals).
 ,w_copy(row_length,rows,0:model_levels) & ! vertical velocity (m/s)
 ,w_max(row_length,rows)                 & ! col max vertical velocity (m/s)
 ,q(row_length,rows,model_levels)        & ! water vapour (kg/kg)
 ,theta(tdims%i_start:tdims%i_end,       & ! Theta (Kelvin)
        tdims%j_start:tdims%j_end,       &
                    1:tdims%k_end)
LOGICAL, INTENT(IN) :: &
  land_mask(row_length, rows)  ! T If land, F ELSEwhere.

REAL, INTENT(IN) ::                         &
  timestep                   & ! timestep (seconds).
 ,frac_land(nunstable)       & ! fraction of land in gridbox
 ,tv1_sd(nunstable)          & ! Approx to standard dev of level of
                               ! 1 virtual temperature (K).
 ,bl_vscale2(nunstable)        ! Velocity scale squared for boundary
                               ! layer eddies (m/s)

! History of convection prognostics - set but not used at present
REAL, INTENT(IN) ::              &
  deep_flag(row_length,rows)     & ! 0-1.0, 1 if deep last time step
 ,past_precip(row_length,rows)   & ! convective precip rate last step
                                   ! or a decayed value.
 ,past_conv_ht(row_length,rows)    ! Convective cloud top on last step
                                   ! (m) (NOT USED at present)

! Additional variables for SCM diagnostics which are dummy in full UM
INTEGER, INTENT(IN) :: &
  nSCMDpkgs          ! No of SCM diagnostics packages

LOGICAL, INTENT(IN) :: &
  L_SCMDiags(nSCMDpkgs) ! Logicals for SCM diagnostics packages


INTEGER, INTENT(INOUT) ::        &
  ntml(row_length,rows)          & ! Number of model levels in the
                                   ! turbulently mixed layer.
 ,ntpar(row_length,rows)         & ! Max levels for parcel ascent
 ,nlcl(row_length,rows)            ! No. of model layers below the
                                   ! lifting condensation level.

LOGICAL, INTENT(INOUT) ::        &
  cumulus(row_length,rows)       & ! Logical indicator for convection
 ,l_shallow(row_length,rows)     & ! Logical indicator for shallow Cu
 ,l_congestus(row_length,rows)   & ! Logical indicator for congestus Cu
 ,l_congestus2(row_length,rows)    ! Logical ind 2 for congestus Cu

REAL, INTENT(INOUT) ::           &
  zh(row_length,rows)              ! Height above surface of top of boundary
                                   ! layer (metres).

REAL, INTENT(OUT) ::             &
  zhpar(row_length,rows)         & ! Height of max parcel ascent (m)
 ,z_lcl(row_length,rows)         & ! Height of lifting condensation
                                   !  level (not a model level) (m)
 ,z_lcl_uv(row_length,rows)      & ! Height of lifting condensation
                                   ! level on nearest uv level (m)
 ,dzh(row_length,rows)           & ! Inversion thickness (m)
 ,qcl_inv_top(row_length,rows)   & ! Parcel water content at inversion top
 ,delthvu(row_length,rows)       & ! Integral of undilute parcel buoyancy
                                   ! over convective cloud layer
                                   ! (for convection scheme)
 ,ql_ad(row_length,rows)         & ! adiabatic liquid water content at
                                   ! inversion or cloud top (kg/kg)
 ,cape(row_length, rows)         & ! CAPE from undilute parcel ascent (m2/s2)
 ,cin(row_length, rows)          & ! CIN from undilute parcel ascent (m2/s2)
 ,entrain_coef(row_length,rows)  & ! Entrainment coefficient
 ,qsat_lcl(row_length,rows)        ! qsat at cloud base


!-----------------------------------------------------------------------
! Local variables
!-----------------------------------------------------------------------
INTEGER ::   &
  i,j        & ! Local Loop counter (horizontal field index).
 ,ii         & ! Local compressed array counter.
 ,k          & ! Local Loop counter (vertical level index).
 ,ntpar_l    & ! cloud top level
 ,mbl          ! Maximum number of model layers allowed in the
               ! mixing layer; set to bl_levels-1.

INTEGER ::   &
  error        ! error code

! Arrays holding various key model levels

INTEGER ::   &
  kshmin(nunstable)        & ! Position of buoyancy minimum above
                             ! topbl (used for shallow Cu diag)
 ,kcucheck(nunstable)      & ! Position of level just below 2.5km
                             ! (used for gradient check to diagnose convection)
 ,k_plume(nunstable)       & ! start level for surface-driven plume
 ,k_max(nunstable)         & ! level of max parcel buoyancy
 ,k_max_dil(nunstable)     & ! level of max parcel buoyancy - dilute
 ,nlcl_min(nunstable)      & ! minimum allowed level of LCL
 ,k_neutral(nunstable)     & ! level of neutral parcel buoyancy
 ,k_neutral_dil(nunstable) & ! level of neutral parcel buoyancy - dilute
 ,k_inv(nunstable)         & ! level from inversion testing
 ,freeze_lev(nunstable)      ! freezing level


! Uncompressed arrays - all points

REAL ::                      &
  z_lcl_nlcl(row_length,rows)   ! Height of nlcl model level (m)

! Compressed arrays store only values for unstable columns
! names as original array plus _c

INTEGER ::                &
  ntml_c(nunstable)       &
 ,nlcl_c(nunstable)       &
 ,nlfc_c(nunstable)       &
 ,ntml_save(nunstable)    &
 ,ntinv_c(nunstable)

REAL ::                                                      &
  q_c(nunstable, model_levels)                               &
 ,qcl_c(nunstable, model_levels)                             &
 ,qcf_c(nunstable, model_levels)                             &
 ,z_full_c(nunstable, model_levels)                          &
 ,z_half_c(nunstable, model_levels)                          &
 ,exner_theta_levels_c(nunstable, model_levels)              &
 ,exner_rho_c(nunstable, model_levels)                       &
 ,P_theta_lev_c(nunstable, model_levels)                     &
 ,P_c(nunstable, model_levels)                               &
 ,cloud_fraction_c(nunstable,model_levels)                   &
 ,zh_c(nunstable)                                            &
 ,zh_save(nunstable)                                         &
 ,zh_itop_c(nunstable)                                       &
 ,delthvu_c(nunstable)                                       &
 ,cape_c(nunstable)                                          &
 ,cin_c(nunstable)                                           &
 ,ql_ad_c(nunstable)                                         &
 ,pstar_c(nunstable)                                         &
 ,qsat_lcl_c(nunstable)

REAL ::                  &
  z_lcl_c(nunstable)     & ! LCL height rounded to nearest model level
 ,zlcl_c(nunstable)        ! Different variable exact height of LCL

LOGICAL ::               &
  l_dilute                 ! if true also DO a dilute ascent


! Arrays only used for unstable calculations

LOGICAL ::                &
  shmin(nunstable)        & ! Flag for finding min in parcel buoyancy below
                            ! 3km (for shallow Cu)
 ,cumulus_c(nunstable)    & ! compressed cumulus
 ,shallow_c(nunstable)      ! compressed shallow

REAL ::                                 &
  t(nunstable, model_levels)            & ! temperature (from theta)
 ,tl(nunstable, model_levels)           & ! Ice/liquid water temperature,
                                          ! but replaced by T in LS_CLD.
 ,qw(nunstable, model_levels)           & ! Total water content
 ,svl(nunstable, model_levels)          & ! Liquid/frozen water virtual
                                          ! static energy over cp.
 ,env_svl(nunstable,model_levels)       & ! Density (virtual) static energy
                                          ! over cp for layer.
 ,par_svl(nunstable,model_levels)       & ! Density (virtual) static energy
                                          ! over cp of parcel for level.
 ,par_svl_dil(nunstable,model_levels)     ! Density (virtual) static energy
                                          ! over cp of parcel for level.

REAL ::                                 &
  t_lcl(nunstable)             & ! Temperature at lifting condensation level.
 ,p_lcl(nunstable)             & ! Pressure at lifting condensation level.
 ,sl_plume(nunstable)          & ! Liquid/frozen water static energy
                                 ! over cp for a plume rising without
                                 ! dilution from level 1.
 ,qw_plume(nunstable)          & ! qw for a plume rising without
                                 ! dilution from level 1.
 ,Dt_dens_parc_T(nunstable)    & ! t_dens_parc-t_dens_env at ntpar
 ,Dt_dens_parc_Tmin(nunstable) & ! t_dens_parc-t_dens_env at kshmin
 ,thv_pert(nunstable)          & ! threshold thv of parcel
 ,delthvu_add(nunstable)       & ! adjustment applied for inversion top
! Added for improved parcel top - mainly used for finding an ascent
! capped by an inversion.
 ,Dt_dens_parc_T2(nunstable)   & ! 2nd copy of Dt_dens_parc_T
 ,delthvu2(nunstable)          & !  2nd copy
 ,zh2(nunstable)               & !  2nd copy
 ,max_buoy(nunstable)          & ! max parcel buoyancy
 ,max_buoy_dil(nunstable)      & ! max parcel buoyancy dilute parcel
 ,ql_ad2(nunstable)            & ! ql_ad 2nd copy
 ,pot_en(nunstable)              ! parcel potential energy when overshooting

! parcel calculation

REAL ::                                &
  T_parc(nunstable, model_levels)      & ! Temperature of parcel.
 ,ql_parc_c(nunstable,model_levels)    & ! parcel water
 ,t_dens_env(nunstable, model_levels)  & ! Density potential temperature
                                         ! of environment.
 ,denv_bydz(nunstable, model_levels)   & ! Gradient of density potential
                                         ! temp in the environment.
 ,dpar_bydz(nunstable, model_levels)   & ! Gradient of density potential
                                         ! temperature of the parcel.
 ,buoy_use(nunstable, model_levels)    & ! parcel buoyancy actually used (K)
 ,buoyancy(nunstable, model_levels)    & ! undilute parcel buoyancy (K)
 ,dqsatdz(nunstable, model_levels)       ! dqsat/dz along adiabat

! Arrays added for dilute parcel calculation

REAL ::                                    &
  entrain_fraction(nunstable,model_levels) & ! fraction of environmental
                                             ! air to mix with parcel
 ,t_parc_dil(nunstable,model_levels)       & ! dilute parcel temeperature
 ,buoyancy_dil(nunstable,model_levels)       ! dilute parcel buoyancy (K)

REAL ::           &
  z_surf          &  ! approx height of top of surface layer
 ,z_neut          &  ! interpolated height of neutral buoyancy
 ,dpe             &  ! change in parcel PE between levels
 ,buoy_ave        &  ! average buoyancy between levels
 ,inc             &  ! CIN/CAPE increment for layer
 ,dz              &  ! Layer depth
 ,a_plume_use     &  ! value of a_plume actually used
 ,factor          &  ! multiplying factor
 ,ent_rate        &  ! entrainment rate
 ,min_precip      &  ! min precip rate
 ,max_precip      &  ! max precip rate
 ,r_over_a           ! radius of a level over radius of earth

! required for average w calculation

REAL ::                                    &
  dmass_theta(nunstable,model_levels) & ! r**2rho*dr on theta levels
 ,w_avg(nunstable)                    & ! mean w over layer (m/s)
 ,w_avg2(nunstable)                     ! mean w over layer (m/s)

REAL ::                         &
  pstar_w_cape_limit(nunstable)   ! scaled critical vertical velocity

CHARACTER (LEN=*), PARAMETER ::  RoutineName = 'CONV_DIAG_COMP_5A'
CHARACTER (LEN=errormessagelength) :: cmessage        ! error message

!-----------------------------------------------------------------------
! Required for SCM diagnostics
!-----------------------------------------------------------------------
REAL :: TmpScm2d(ScmRowLen,ScmRow)              ! Single level field
REAL :: TmpScm3d(ScmRowLen,ScmRow,model_levels) ! Full field



!
! Model constants:
!

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!-----------------------------------------------------------------------
! Mixing ratio, r,  versus specific humidity, q
!
! In most cases the expression to first order are the same
!
!  tl = T - (lc/cp)qcl - [(lc+lf)/cp]qcf
!  tl = T - (lc/cp)rcl - [(lc+lf)/cp]rcf  - equally correct definition
!
! thetav = theta(1+cvq)         accurate
!        = theta(1+r/repsilon)/(1+r) ~ theta(1+cvr) approximate
!
! svl = (tl+gz/cp)*(1+(1/repsilon-1)qt)
!     ~ (tl+gz/cp)*(1+(1/repsilon-1)rt)
!
! dqsat/dT = repsilon*Lc*qsat/(R*T*T)
! drsat/dT = repsilon*Lc*rsat/(R*T*T)  equally approximate
!
! Only altering the expression for vapour pressure
!
!  e = qp/repsilon       - approximation
!  e = rp/(repsilon+r)   - accurate
!
!-----------------------------------------------------------------------
! 1.0 set variables
!-----------------------------------------------------------------------
!  Set MBL, "maximum number of boundary levels" for the purposes of
!  boundary layer height calculation.

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

mbl = bl_levels - 1

!$OMP  PARALLEL DEFAULT(SHARED)                                         &
!$OMP& PRIVATE(ii, i, j, r_over_a, k, z_surf, a_plume_use)

!-----------------------------------------------------------------------
! 1.1 initialise just unstable points
!-----------------------------------------------------------------------

!$OMP DO SCHEDULE(STATIC)
DO ii=1,nunstable

  cumulus_c(ii)    = .FALSE.
  shallow_c(ii)    = .FALSE.

  kcucheck(ii)   = 1
  freeze_lev(ii) = 1

  ntml_c(ii) = 1
  nlcl_c(ii) = 1
  nlfc_c(ii) = 0  ! implies unset

END DO
!$OMP END DO

!-----------------------------------------------------------------------
! 2.0 Calculate mass of layers for use later
!-----------------------------------------------------------------------

!$OMP DO SCHEDULE(STATIC)
DO k=1, model_levels-1
  DO ii=1, nunstable
    i = index_i(ii)
    j = index_j(ii)
    r_over_a = r_theta_levels(i,j,k)/planet_radius
    dmass_theta(ii,k) = rho_theta(i,j,k)*r_over_a*r_over_a              &
                           *(r_rho_levels(i,j,k+1)-r_rho_levels(i,j,k))
  END DO
END DO
!$OMP END DO

!-----------------------------------------------------------------------
! 3.0 Calculate various quantities required by the parcel calculations
!-----------------------------------------------------------------------
! compress boundary layer depth

!$OMP DO SCHEDULE(STATIC)
DO ii=1,nunstable
  i = index_i(ii)
  j = index_j(ii)
  zh_c(ii)  = zh(i,j)
  z_lcl_c(ii) = z_half(i,j,nlcl_c(ii)+1)
  pstar_c(ii) = pstar(i,j)
END DO
!$OMP END DO

!$OMP DO SCHEDULE(STATIC)
DO k = 1, model_levels
  DO ii=1, nunstable
    i = index_i(ii)
    j = index_j(ii)
    t(ii,k) = theta(i,j,k) * exner_theta_levels(i,j,k)

    ! initialise t_parc at all points
    ! added for safety of qsat_mix calls later

    t_parc(ii,k) = t(ii,k)
    q_c(ii,k)   = q(i,j,k)
    qcl_c(ii,k) = qcl(i,j,k)
    qcf_c(ii,k) = qcf(i,j,k)
    z_full_c(ii,k) = z_full(i,j,k)
    z_half_c(ii,k) = z_half(i,j,k)
    exner_theta_levels_c(ii,k) = exner_theta_levels(i,j,k)
    exner_rho_c(ii,k)          = exner_rho(i,j,k)
    P_c(ii,k)              = p(i,j,k)
    P_theta_lev_c(ii,k)    = P_theta_lev(i,j,k)
    cloud_fraction_c(ii,k) = cloud_fraction(i,j,k)
  END DO
END DO
!$OMP END DO

!-----------------------------------------------------------------------
! 3.1 Calculate total water content, qw and Liquid water temperature, tl
!     Definitions for qw and tl the same whether q, qcl, qcf  are
!     specific humidities  or mixing ratio.
!
!-----------------------------------------------------------------------

!$OMP DO SCHEDULE(STATIC)
DO k=1,model_levels
  DO ii=1, nunstable

    ! Total water   - BL doc P243.10

    qw(ii,k) = q_c(ii,k) + qcl_c(ii,k) + qcf_c(ii,k)

    ! Liquid water temperature as defined  BL doc P243.9

    tl(ii,k) = t(ii,k) - lcrcp*qcl_c(ii,k) - lsrcp*qcf_c(ii,k)


    ! Calculate svl: conserved variable  a form of moist static energy /cp
    !       svl = (tl+gz/cp)*(1+(1/repsilon-1)qt)  - specific humidity

    svl(ii,k) = ( tl(ii,k) + gamma_dry * z_full_c(ii,k) )                   &
                                     * ( 1.0 + c_virtual*qw(ii,k) )

    ! Density potential temperature of environment (K)

    t_dens_env(ii,k)=t(ii,k)*(1.0+c_virtual*q_c(ii,k)-qcl_c(ii,k)-qcf_c(ii,k))

  END DO
END DO
!$OMP END DO

!-----------------------------------------------------------------------
! 4.0 Parts of parcel calculation common to all options
!-----------------------------------------------------------------------
! 4.1 Work out initial parcel properties and LCL
!-----------------------------------------------------------------------

!$OMP DO SCHEDULE(STATIC)
DO ii=1, nunstable

  k_plume(ii) = 1
  !-----------------------------------------------------------------------
  ! Only perform parcel ascent if unstable
  ! Start plume ascent from grid-level above top of surface layer, taken
  ! to be at a height, z_surf, given by 0.1*zh
  !-----------------------------------------------------------------------
  z_surf = 0.1 * zh_c(ii)

  DO WHILE( z_full_c(ii,k_plume(ii))  <   z_surf .AND.                    &
    !                   ! not reached z_surf
                svl(ii,k_plume(ii)+1)  <   svl(ii,k_plume(ii)) )
    !                   ! not reached inversion

    k_plume(ii) = k_plume(ii) + 1

  END DO
END DO          ! loop over ii
!$OMP END DO

a_plume_use = a_plume
IF (forced_cu >= on) a_plume_use = 0.0  ! no need to add bonus buoyancy

IF (limit_pert_opt == 2) THEN

!$OMP DO SCHEDULE(STATIC)
  DO ii=1, nunstable
    sl_plume(ii) = tl(ii,k_plume(ii)) + gamma_dry * z_full_c(ii,k_plume(ii))
    thv_pert(ii) = MIN( MAX( a_plume_use,                                 &
                       MIN( max_t_grad*zh_c(ii), b_plume*tv1_sd(ii) ) ),  &
                       max_diag_thpert )
    qw_plume(ii) = qw(ii,k_plume(ii))
  END DO
!$OMP END DO

ELSE IF (limit_pert_opt == 0 .OR. limit_pert_opt == 1) THEN

!$OMP DO SCHEDULE(STATIC)
  DO ii=1, nunstable
    sl_plume(ii) = tl(ii,k_plume(ii)) + gamma_dry * z_full_c(ii,k_plume(ii))
    thv_pert(ii) = MAX( a_plume_use,                                      &
                       MIN( max_t_grad*zh_c(ii), b_plume*tv1_sd(ii) ) )
    qw_plume(ii) = qw(ii,k_plume(ii))
  END DO
!$OMP END DO

END IF

!$OMP END PARALLEL

!-----------------------------------------------------------------------
! 4.2 Calculate temperature and pressure of lifting condensation level
!       using approximations from Bolton (1980)
!-----------------------------------------------------------------------
!
CALL lift_cond_lev ( nunstable, model_levels, k_plume,                    &
                     pstar_c, q_c, t,                                     &
                     P_theta_lev_c, exner_rho_c, z_half_c,                &
                     t_lcl, p_lcl, zlcl_c, qsat_lcl_c )

!-----------------------------------------------------------------------
! Convection scheme requires NLCL to be at least 2.
! For model stability over mountains, also require ZLCL > 150m
! (approx nlcl=2 for G3 levels)
!-----------------------------------------------------------------------

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ii, k, i, j)

!$OMP DO SCHEDULE(STATIC)
DO ii=1, nunstable

  k=3
  DO WHILE ( z_half_c(ii,k) < 150.0 .AND. k < mbl )
    k=k+1
  END DO
  nlcl_min(ii) = k-1

END DO
!$OMP END DO

!-----------------------------------------------------------------------
! Reset zh  (at this point in the code ntml is initialised as =1)
!-----------------------------------------------------------------------

!$OMP DO SCHEDULE(STATIC)
DO j=1,rows
  DO i=1,row_length
    zh(i,j) = z_half(i,j,ntml(i,j)+1)
  END DO
END DO
!$OMP END DO

!$OMP DO SCHEDULE(STATIC)
DO ii=1, nunstable
  i = index_i(ii)
  j = index_j(ii)
  zh_c(ii) = zh(i,j)
  zh2(ii)  = zh_c(ii)
  z_lcl(i,j)    = zlcl_c(ii)         ! expand up accurate z_lcl
  qsat_lcl(i,j) = qsat_lcl_c(ii)     ! expand up qsat at cloud base
END DO
!$OMP END DO

!-----------------------------------------------------------------------
! Find NLCL
!-----------------------------------------------------------------------
!
!    ---------------   p,T,q       nlcl+1  , p_theta(nlcl+2)
!
!    - - - - - - - -   uv,p,rho    nlcl+1,  z_lcl , p(nlcl+1)     either
!     + + + + + + + +   lcl, Plcl, not a model level            lower part
!    ---------------   p,T,q       nlcl , p_theta(nlcl+1)          of layer
!
!    - - - - - - - -   uv,p,rho    nlcl   p(nlcl)
!
!-----------------------------------------------------------------------
!
!    ---------------   p,T,q       nlcl+1  , p_theta(nlcl+2)
!     + + + + + + + +   lcl, Plcl, not a model level
!
!    - - - - - - - -   uv,p,rho    nlcl+1,  z_lcl , p(nlcl+1)        or
!                                                                upper part
!    ---------------   p,T,q       nlcl , p_theta_lev(nlcl+1)      of layer
!
!    - - - - - - - -   uv,p,rho    nlcl   p(nlcl)
!
!-----------------------------------------------------------------------


DO k = 2,model_levels
!$OMP DO SCHEDULE(STATIC)
  DO ii=1, nunstable

    ! NLCL level
    IF ( p_lcl(ii)  <   P_c(ii,k) ) THEN

      ! compressed copies
      nlcl_c(ii)  = k-1
      z_lcl_c(ii) = z_half_c(ii,nlcl_c(ii)+1)
    END IF     ! test on p_lcl

    ! Find level just below 2.5km (for use in Cu diagnosis)

    IF ( z_full_c(ii,k)  >   2500.0                                      &
               .AND. kcucheck(ii)  ==  1 ) kcucheck(ii) = k-1

    ! Freezing level

    IF (t(ii,k) <  tm .AND. t(ii,k-1) >= tm) THEN
      IF (freeze_lev(ii) == 1) THEN
        freeze_lev(ii) = k
      END IF
    END IF

  END DO       ! ii loop
!$OMP END DO
END DO         ! k loop

! expand to full arrays
!$OMP DO SCHEDULE(STATIC)
DO ii=1, nunstable
  i = index_i(ii)
  j = index_j(ii)
  nlcl(i,j) = nlcl_c(ii)
  z_lcl_nlcl(i,j) = z_half_c(ii,nlcl_c(ii)+1)
  z_lcl_uv(i,j)   = z_full_c(ii,nlcl_c(ii))
END DO       ! ii loop
!$OMP END DO

!$OMP END PARALLEL

!=======================================================================
! 5.0 Parcel calculation - various options available
!                          (only on unstable points)
!=======================================================================
!
!  icvdiag option   |      explaination
!  -------------------------------------------------------------------
!     1             |  undilute parcel calculation
!  -  > 2 - - - -   |  - dilute parcel calculation - - - - - - - - -
!     2             | 0.55/z entrainment rates
!     3             |  1/z entrainment rate
!  -------------------------------------------------------------------
!     >3        Experimental options
!     4             |  Diurnal cycle over land entrainment rates,
!                   !  Ocean used 0.55/Z dilution of parcel
!     5             |  Diurnal cycle over land entrainment rates,
!                   !  Ocean undilute ascent of parcel
!=======================================================================
! Original 5A code - undilute parcel
!=======================================================================
IF (icvdiag == 1) THEN
  !=======================================================================

  l_dilute = .FALSE.

  !=======================================================================
  ! Dilute parcel ascent options
  !=======================================================================
ELSE IF (icvdiag >= 2 .AND. icvdiag < 6 ) THEN

  l_dilute = .TRUE.       ! dilute ascent required as well as undilute

  ! Set top wet level entrainment to zero
  k=model_levels
  DO ii=1, nunstable
    entrain_fraction(ii,k) = 0.0
  END DO

  !-----------------------------------------------------------------------
  !  Option 2 - 0.55/z entrainment rate as Jakob & Siebesma
  !             Entrainment rate not passed to convection scheme.
  !-----------------------------------------------------------------------

  IF (icvdiag == 2) THEN

    DO k= 1,model_levels-1
      DO ii=1, nunstable
        i = index_i(ii)
        j = index_j(ii)
        dz = r_rho_levels(i,j,k+1)-r_rho_levels(i,j,k)
        entrain_fraction(ii,k) = 0.55 * dz/z_full_c(ii,k)
      END DO
    END DO

    !-----------------------------------------------------------------------
    !  Option 3 - 1/z entrainment rate, similar to many CRM results for
    !             deep and shallow convection.
    !             Entrainment rate not passed to convection scheme.
    !-----------------------------------------------------------------------
  ELSE IF (icvdiag == 3) THEN

    DO k= 1,model_levels-1
      DO ii=1, nunstable
        i = index_i(ii)
        j = index_j(ii)
        dz = r_rho_levels(i,j,k+1)-r_rho_levels(i,j,k)
        entrain_fraction(ii,k) = 1.0 * dz/z_full_c(ii,k)
      END DO
    END DO

    !-----------------------------------------------------------------------
    !  Option 4 - 0.55/z entrainment rate  for ocean
    !             Diurnal varying rate over land.
    !             Entrainment rate for land passed to convection scheme.
    !-----------------------------------------------------------------------

  ELSE IF (icvdiag == 4 ) THEN

    !--------------------------------------------------------------------
    ! Diurnal cycle changes
    !--------------------------------------------------------------------
    ! Only applied to land points (points with > 80% land if coastal).
    ! x/z entrainment rates - where x depends on BL depth (zlcl)
    ! If zlcl > 1000m THEN assumes convection uses so called equilbrium
    ! entrainment rates for land, i.e. 0.55/z ?
    ! No constraint on maximum entrainment rate (formula => 10.23,
    ! if zlcl=0.0)
    !--------------------------------------------------------------------

    DO k= 1,model_levels-1
      DO ii=1, nunstable
        i = index_i(ii)
        j = index_j(ii)

        IF (frac_land(ii) > 0.8) THEN
          ent_rate = 0.55 +                                                 &
                        8.0*(1.2-(z_lcl(i,j)/1000.0))*(1.2-(z_lcl(i,j)/1000.0))
          entrain_coef(i,j) = ent_rate

          ! While very high entrainment rates are ok in the diagnosis
          ! the shallow convection scheme does not perform well
          ! if passed an entrainment rate ~10. Therefore restricting
          ! entrainment rates passed to convection.
          IF (entrain_coef(i,j) > 3.5) THEN
            entrain_coef(i,j) = 3.5
          END IF
          IF (z_lcl(i,j) > 1200.0) THEN
            ent_rate = 0.55
            entrain_coef(i,j) = 0.55
          END IF

          !--------------------------------------------------------------------
          ! Suppressing convection over land where there is steep orography
          ! and very high vertical velocities may tend to cause gridpoint
          ! storms. A test has been added to revert to the standard sea
          ! entrainment rates where the w based CAPE closure is being used to
          ! control gridpoint storms.
          !--------------------------------------------------------------------

          IF (cldbase_opt_dp == 3 .OR. cldbase_opt_dp == 4 .OR.             &
              cldbase_opt_dp == 5 .OR. cldbase_opt_dp == 6 .OR.             &
              cldbase_opt_dp == 7 .OR. cldbase_opt_dp == 8 ) THEN
            ! Pressure weighted w_max test
            pstar_w_cape_limit(ii) = w_cape_limit * pstar_c(ii) / 1.0e5
            IF (w_max(i,j) > pstar_w_cape_limit(ii)) THEN
              ent_rate = 0.55
              entrain_coef(i,j) = -99.0
            END IF
          END IF

        ELSE               ! Sea points  or coastal points
          !------------------------------------------------------------------
          ! Sea points - using a dilute 0.55/z diagnosis but standard
          !              entrainment rates for shallow and deep convection
          !------------------------------------------------------------------
          ent_rate = 0.55
          entrain_coef(i,j) = -99.0

        END IF        ! land fraction test

        dz = r_rho_levels(i,j,k+1)-r_rho_levels(i,j,k)
        entrain_fraction(ii,k) = ent_rate * dz/z_full_c(ii,k)

      END DO          ! unstable point loop
    END DO            ! level loop

    !-----------------------------------------------------------------------
    !  Option 5 - undilute ocean diagnosis
    !             Diurnal varying dilute ascent over land.
    !             Entrainment rate for land passed to convection scheme.
    !-----------------------------------------------------------------------

  ELSE IF (icvdiag == 5 ) THEN

    !--------------------------------------------------------------------
    ! Diurnal cycle changes
    !--------------------------------------------------------------------
    ! Only applied to land points (points with > 80% land if coastal).
    ! x/z entrainment rates - where x depends on BL depth (zlcl)
    ! If zlcl > 1000m THEN assumes convection uses so called equilbrium
    ! entrainment rates for land, i.e. 0.55/z ?
    ! No constraint on maximum entrainment rate (formula => 10.23,
    ! if zlcl=0.0)
    !--------------------------------------------------------------------

    DO k= 1,model_levels-1
      DO ii=1, nunstable
        i = index_i(ii)
        j = index_j(ii)

        IF (frac_land(ii) > 0.8) THEN
          ent_rate = 0.55 +                                                 &
                        8.0*(1.2-(z_lcl(i,j)/1000.0))*(1.2-(z_lcl(i,j)/1000.0))
          entrain_coef(i,j) = ent_rate

          ! While very high entrainment rates are ok in the diagnosis
          ! the convection scheme does not perform well
          ! if passed an entrainment rate ~10. Therefore restricting
          ! entrainment rates passed to convection.
          IF (entrain_coef(i,j) > 3.5) THEN
            entrain_coef(i,j) = 3.5
          END IF
          IF (z_lcl(i,j) > 1200.0) THEN
            ent_rate = 0.55
            entrain_coef(i,j) = 0.55   ! uses 0.55/z
          END IF

          !--------------------------------------------------------------------
          ! Suppressing convection over land where there is steep orography
          ! and very higher vertical velocities may tend to cause gridpoint
          ! storms. A test has been added to revert to the standard sea
          ! entrainment rates where the w based CAPE closure is being used to
          ! control gridpoint storms.
          !--------------------------------------------------------------------

          IF (cldbase_opt_dp == 3 .OR. cldbase_opt_dp == 4 .OR.             &
              cldbase_opt_dp == 5 .OR. cldbase_opt_dp == 6 .OR.             &
              cldbase_opt_dp == 7 .OR. cldbase_opt_dp == 8 ) THEN
            ! Pressure weighted w_max test
            pstar_w_cape_limit(ii) = w_cape_limit * pstar_c(ii) / 1.0e5
            IF (w_max(i,j) > pstar_w_cape_limit(ii)) THEN
              ent_rate = 0.55
              entrain_coef(i,j) = -99.0
            END IF
          END IF

        ELSE               ! Sea points  or coastal points
          !------------------------------------------------------------------
          ! Sea points - no dilute ascent
          !------------------------------------------------------------------
          ent_rate = 0.0
          entrain_coef(i,j) = -99.0

        END IF        ! land fraction test

        dz = r_rho_levels(i,j,k+1)-r_rho_levels(i,j,k)
        entrain_fraction(ii,k) = ent_rate * dz/z_full_c(ii,k)

      END DO          ! unstable point loop
    END DO            ! level loop

  END IF    !   diagnosis option

ELSE

  error = 2     ! fatal return code
  ! Statement assumes icvdiag is in the range that format statement covers.
  WRITE(cmessage,'(a40,i6)')                                   &
     'Convection diagnosis option not allowed ',icvdiag
  CALL ereport(routinename,error,cmessage)

  !=======================================================================
  !    End of choice of diagnosis code
  !=======================================================================
END IF

!-----------------------------------------------------------------------
! Parcel ascent
!-----------------------------------------------------------------------
! Parcel starts from level k_plume and is lifted up.
! If a dilute parcel ascent - mix in environmental air above lifting
! condensation level
!----------------------------------------------------------------------
CALL parcel_ascent ( nunstable, nSCMDpkgs,                             &
                       nlcl_c, k_plume,                                &
                       l_dilute, L_SCMDiags,                           &
                       sl_plume, qw_plume,                             &
                       t, tl, qcl_c, qcf_c, qw, t_dens_env,            &
                       p_theta_lev_c, exner_theta_levels_c,            &
                       z_full_c, entrain_fraction,                     &
                       T_parc, t_parc_dil, ql_parc_c,                  &
                       buoyancy,buoyancy_dil,                          &
                       env_svl, par_svl, par_svl_dil,                  &
                       denv_bydz, dpar_bydz, dqsatdz )


!-----------------------------------------------------------------------
! SCM diagnostics if no unstable points
!-----------------------------------------------------------------------
IF (model_type == mt_single_column) THEN
  IF (nunstable == 0 .AND. l_dilute) THEN

    !  zero output
    TmpScm3d(:,:,:) = 0.0

    IF (l_scmdiags(scmdiag_conv)) THEN
      CALL scmoutput(TmpScm3d,'ql_parc',                               &
           'parcel water                ','kg/kg',                     &
           t_inst,d_wet,default_streams,'',routinename)
      CALL scmoutput(TmpScm3d,'ql_parc_dil',                           &
           'parcel water dilute plume','kg/kg',                        &
           t_inst,d_wet,default_streams,'',routinename)
      CALL scmoutput(TmpScm3d,'sl_parc',                               &
           'Parcel energy          ','K',                              &
           t_inst,d_wet,default_streams,'',routinename)
      CALL scmoutput(TmpScm3d,'qw_parc',                               &
           'total parcel water          ','kg/kg',                     &
           t_inst,d_wet,default_streams,'',routinename)
    END IF ! scmdiag_conv

  END IF ! test on no unstable points

END IF ! model_type

!-----------------------------------------------------------------------
! Tests on parcel ascent - look for an inversion or not
!-----------------------------------------------------------------------

SELECT CASE(cvdiag_inv)

CASE (0)    ! no inversion testing

  !-----------------------------------------------------------------------
  !   Now compare plume s_VL with each model layer s_VL in turn to
  !     find the first time that plume has negative buoyancy.
  !-----------------------------------------------------------------------
  CALL cv_parcel_neutral_dil(nunstable,                                     &
             nlcl_c,k_plume,l_dilute,                                       &
             z_lcl_c, thv_pert, z_full_c, z_half_c, exner_theta_levels_c,   &
             buoyancy, buoyancy_dil, t_dens_env, dqsatdz,                   &
             zh_c,                                                          &
             k_max,k_max_dil,k_neutral,k_neutral_dil,                       &
             max_buoy,max_buoy_dil,ql_ad_c,delthvu_c,                       &
             cape_c,cin_c)

  !-----------------------------------------------------------------------
  ! Default parcel top properties are assumed to be those when the
  ! ascent reaches the level of neutral buoyancy.
  !-----------------------------------------------------------------------

  IF (l_dilute) THEN

    DO ii=1,nunstable
      ntml_c(ii)    = k_neutral_dil(ii)
      ntinv_c(ii)   = ntml_c(ii)
      zh_itop_c(ii) = -1.0   ! to test when set sensibly
    END DO
    DO k=1,model_levels
      DO ii=1,nunstable
        buoy_use(ii,k)=buoyancy_dil(ii,k)+thv_pert(ii)
      END DO
    END DO
  ELSE     ! undilute ascent

    DO ii=1,nunstable
      ntml_c(ii)    = k_neutral(ii)
      ntinv_c(ii)   = ntml_c(ii)
      zh_itop_c(ii) = -1.0   ! to test when set sensibly
    END DO    ! ii loop
    DO k=1,model_levels
      DO ii=1,nunstable
        buoy_use(ii,k)=buoyancy(ii,k)+thv_pert(ii)
      END DO
    END DO    ! ii loop

  END IF  ! test on l_dilute

  !-----------------------------------------------------------------------
  ! Estimate inversion thickness
  !-----------------------------------------------------------------------
  IF (forced_cu >= on .OR. bl_res_inv == on) THEN

    ! find the level of free convection
    !  - nlfc is the first grid level above the LCL that is buoyant
    DO k=2,model_levels-1
      DO ii=1,nunstable
        IF (nlfc_c(ii) == 0 .AND. k > nlcl_c(ii) .AND.                    &
            ! LFC unset but above LCL
            buoy_use(ii,k) > 0.0) THEN
            ! become positively buoyant
          nlfc_c(ii) = k
        END IF
      END DO
    END DO

    DO ii=1,nunstable
      delthvu_add(ii)=0.0   ! initiailise
      k = ntml_c(ii)
      ! Find height where buoy_use=0, linearly interpolating
      z_neut = z_full_c(ii,k) - buoy_use(ii,k) *                          &
                        (z_full_c(ii,k+1)-z_full_c(ii,k))/                &
                        (buoy_use(ii,k+1)-buoy_use(ii,k))
      dz = z_full_c(ii,k+1) - z_neut
      ! Note buoy_use=0 at z_neut, hence integrating to z_full(k+1) gives:
      pot_en(ii) = -0.5*buoy_use(ii,k+1)*dz*g/t_dens_env(ii,k+1)
      IF (pot_en(ii) > bl_vscale2(ii) ) THEN
        zh_itop_c(ii) = MAX( zh_c(ii), z_neut -                           &
            2.0*bl_vscale2(ii)*t_dens_env(ii,k+1)/(buoy_use(ii,k+1)*g) )
      END IF
    END DO

    DO k=2,model_levels-1
      DO ii=1,nunstable
        IF (zh_itop_c(ii) < 0.0 .AND. k >= ntml_c(ii)+1) THEN
          ! look for inversion top
          dz = z_full_c(ii,k+1) - z_full_c(ii,k)
          buoy_ave = 0.5*(buoy_use(ii,k)+buoy_use(ii,k+1))
          dpe = - buoy_ave * dz * g / t_dens_env(ii,k)
          IF ( pot_en(ii)+dpe > bl_vscale2(ii) ) THEN
            zh_itop_c(ii) = MAX( zh_c(ii), z_full_c(ii,k) +               &
                         (pot_en(ii)-bl_vscale2(ii))*t_dens_env(ii,k)/    &
                         ( buoy_ave*g ) )
            IF (zh_itop_c(ii) >= z_half_c(ii,k+1)) THEN
              ntinv_c(ii) = k  ! marks top level within inversion
            ELSE
              ntinv_c(ii) = k-1  ! marks top level within inversion
            END IF
            IF (delthvu_add(ii) > 0.0) THEN
              ! if parcel has overcome inhibition to become positively
              ! buoyant then adjust delthvu to allow convection triggering
              delthvu_c(ii)=delthvu_c(ii)+delthvu_add(ii)
            END IF
          ELSE
            pot_en(ii) = pot_en(ii) + dpe
            dz = z_half_c(ii,k+1) - z_half_c(ii,k)
            delthvu_add(ii) = delthvu_add(ii) + buoy_use(ii,k)*dz         &
                                          /exner_theta_levels_c(ii,k)
          END IF
        END IF
      END DO  ! over ii
    END DO  ! over k

    DO ii=1,nunstable
      ntml_save(ii) = ntml_c(ii)
      zh_save(ii)   = zh_c(ii)
      IF (ntinv_c(ii) > nlfc_c(ii)) THEN
        ! Inversion top is above LFC so implies free rather than forced cu
        !  - use inversion top as top of parcel ascent in cumulus diagnosis
        ntml_c(ii) = ntinv_c(ii)
        zh_c(ii)   = zh_itop_c(ii)
        ! Test whether cloud within inversion is capped by a significant
        ! inversion (lapse rate greater than moist adiabatic, given
        ! by dpar_bydz) - if not, trigger cumulus
        k = ntinv_c(ii)
        IF ( ntinv_c(ii) > nlcl_c(ii) .AND.                             &
             denv_bydz(ii,k+1) < 1.25*dpar_bydz(ii,k+1) ) THEN
          cumulus_c(ii) = .TRUE.
        END IF
      END IF   ! ntinv above lfc

    END DO

  END IF  ! test on forced_cu >= on or bl_res_inv == on

  !-----------------------------------------------------------------------
  ! Average vertical velocity over a layer  - required for shallow
  !   convection test.
  !-----------------------------------------------------------------------
  ! Layer from top in cloud w to value 1500km above cloud top
  !-----------------------------------------------------------------------

  CALL mean_w_layer(nunstable,row_length,rows,model_levels,                 &
             ntml_c, index_i, index_j,                                      &
             1500.0, z_full_c, z_half_c, w_copy, dmass_theta,               &
             w_avg)


CASE (1)    ! Original inversion test for shallow convection
           ! Only available for undilute ascent

  !-----------------------------------------------------------------------
  !   Now compare plume s_VL with each model layer s_VL in turn to
  !     find the first time that plume has negative buoyancy.
  !-----------------------------------------------------------------------

  CALL cv_parcel_neutral_inv(nunstable,                                     &
             nlcl_c,k_plume,                                                &
             z_lcl_c, thv_pert, z_full_c, z_half_c, exner_theta_levels_c,   &
             buoyancy, t_dens_env, dqsatdz, denv_bydz, dpar_bydz,           &
             zh_c, zh2,                                                     &
             k_max,k_neutral,k_inv, kshmin, shmin,                          &
             max_buoy,dt_dens_parc_t,ql_ad_c,delthvu_c,cape_c,cin_c,        &
             dt_dens_parc_t2,dt_dens_parc_tmin,ql_ad2,                      &
             delthvu2)

  !-----------------------------------------------------------------------
  ! Average vertical velocity over a layer  - required for shallow
  !   convection test.
  !-----------------------------------------------------------------------
  ! Layer from top in cloud w to value 1500km above cloud top
  !-----------------------------------------------------------------------

  CALL mean_w_layer(nunstable,row_length,rows,model_levels,                 &
             k_neutral, index_i, index_j,                                   &
             1500.0, z_full_c, z_half_c, w_copy, dmass_theta,               &
             w_avg)

  CALL mean_w_layer(nunstable,row_length,rows,model_levels,                 &
             k_inv, index_i, index_j,                                       &
             1500.0, z_full_c, z_half_c, w_copy, dmass_theta,               &
             w_avg2)

  !-----------------------------------------------------------------------
  ! Default parcel top properties are assumed to be those when the
  ! ascent reaches the level of neutral buoyancy. These may not be those
  ! required in the case of shallow convection.
  ! Shallow convection requires the possible identification of an inversion
  ! at the top of the ascent. This may not be detected by the LNB test.
  ! The gradient tests are designed to detect the shallow top.
  !-----------------------------------------------------------------------
  ! Modify top if ascent is likely to be shallow

  DO ii=1,nunstable

    IF (shmin(ii) ) THEN    ! found an inversion
      ! points where k_inv not the same as k_neutral and level below freezing
      ! may be shallow or congestus or deep

      IF (k_inv(ii) == k_neutral(ii)) THEN
        !  Both methods give same answer for top level leave shmin set
        ntml_c(ii) = k_neutral(ii)

        ! Inversion top lower than level of neutral buoyancy.
        ! Check also, either below freezing level or less than 2500m for
        ! shallow convection.

      ELSE IF ((k_inv(ii) <  freeze_lev(ii) .OR.                       &
                         z_full_C(ii,k_inv(ii)+1)  <=  2500.0 )        &
                 .AND. k_inv(ii) <  k_neutral(ii) ) THEN


        IF ( (z_full_c(ii,kshmin(ii)) - z_full_c(ii,k_inv(ii)))        &
             <=  1.25*(z_half_c(ii,k_inv(ii)+1) - z_lcl_c(ii)) .AND.    &
             (dt_dens_parc_tmin(ii)  <=  0.55*dt_dens_parc_t2(ii))     &
             .AND.     (w_avg2(ii)  <   0.0)  ) THEN

          ! May be shallow or congestus
          ! set values to those found from inversion testing
          ntml_c(ii)  = k_inv(ii)
          delthvu_c(ii) = delthvu2(ii)
          zh_c(ii)    = zh2(ii)
          w_avg(ii)   = w_avg2(ii)
          ql_ad_c(ii) = ql_ad2(ii)
          dt_dens_parc_t(ii) = dt_dens_parc_t2(ii)

        ELSE   ! Assume not shallow or congestus

          ntml_c(ii) = k_neutral(ii)
          shmin(ii) = .FALSE.  ! inversion top found not good
                                 ! don't DO shallow tests
        END IF

      ELSE   ! Assume deep  and therefore top LNB
        ntml_c(ii) = k_neutral(ii)
        shmin(ii) = .FALSE.  ! inversion top found not good
                                 ! don't DO shallow tests
      END IF        ! tests on k_inv

    ELSE    !  No inversion found  i.e. shmin=false

      ntml_c(ii) = k_neutral(ii)

    END IF     ! shmin test

  END DO    ! ii loop

END SELECT  ! test on type of testing inversion of not.


!-----------------------------------------------------------------------
! Test which points are cumulus
!-----------------------------------------------------------------------
CALL cumulus_test (nunstable,mbl,                                      &
                      ntml_c, k_plume, kcucheck,                       &
                      zh_c, frac_land, qw, cloud_fraction_c, z_full_c, &
                      z_lcl_c, nlcl_c, cumulus_c )

!==============================================================================
! Original shallow/deep diagnosis - No congestus diagnosed
! or a run with no parametrized convection so iconv_congestus == imdi 

IF (iconv_congestus == 0 .OR. iconv_congestus == imdi) THEN

  !==============================================================================

  SELECT CASE(cvdiag_inv)

  CASE (0)    ! no inversion testing
    !CDIR NODEP
    DO ii=1,nunstable
      i = index_i(ii)
      j = index_j(ii)
      IF ( cumulus_c(ii) ) THEN

        !---------------------------------------------------------------------
        ! If cumulus has been diagnosed, determine whether it is shallow
        ! or deep convection
        !---------------------------------------------------------------------
        ! Conditions for shallow convection
        !  (1) top of parcel ascent < 2500. or T (top of parcel) > TM
        !  (2) Additional condition   w_avg < cvdiag_sh_wtest
        !---------------------------------------------------------------------
        ntpar_l = ntml_c(ii)

        IF ( z_full_c(ii,ntpar_l)  <=  2500.0 .OR.                      &
                 t(ii,ntpar_l)  >=  tm ) THEN

          IF (l_dilute .OR. iconv_deep == 2) THEN
            IF (frac_land(ii) > 0.8) THEN      ! land points no w test

              shallow_c(ii) = .TRUE.

              ! If undilute parcel would have been deep then reset to go
              ! through deep scheme but using the new entrain rate.

              IF (t(ii,k_neutral(ii)) < tm) THEN
                shallow_c(ii) = .FALSE.
              END IF

            ELSE                       ! ocean points do w test

              IF (w_avg(ii)  <  cvdiag_sh_wtest) THEN
                shallow_c(ii) = .TRUE.
              END IF
            END IF

          ELSE   ! undilute ascent (icvdiag=1)

            IF ( w_avg(ii)  <  cvdiag_sh_wtest ) THEN
              shallow_c(ii) = .TRUE.
            END IF

          END IF  ! l_dilute

        END IF  !  height and temp test

      END IF        ! test on cumulus

    END DO          ! ii loop


  CASE (1)    ! Orignal inversion test
    !CDIR NODEP
    DO ii=1,nunstable
      i = index_i(ii)
      j = index_j(ii)
      IF ( cumulus_c(ii) ) THEN

        !---------------------------------------------------------------------
        ! If cumulus has been diagnosed, determine whether it is shallow
        ! or deep convection
        !---------------------------------------------------------------------
        ! Conditions for shallow convection
        ! (1) w_avg < cvdiag_sh_wtest     (descending air or weak ascent)
        ! (2) top of parcel ascent < 2500. or T (top of parcel) > TM
        ! (3) height of min buoyancy (above Bl) - height of parcel top T level
        !      <1.25(height parcel top - z lifting condensation level)
        ! (4) t_dens_parc -t_dens_env at kshmin <0.55t_dens_parc -t_dens_env
        !     at ntpar
        !
        ! The last 2 conditions are looking for a strong inversion at the top
        ! of the shallow cumulus.
        !---------------------------------------------------------------------
        IF ( shmin(ii) ) THEN
          ntpar_l = ntml_c(ii)
          ! New deep turbulence scheme no w test over land for shallow
          IF (iconv_deep == 2 .AND. frac_land(ii) > 0.8) THEN
            IF ( (z_full_c(ii,ntpar_l)  <=  2500.0 .OR.                   &
                     t(ii,ntpar_l)  >=  tm)                               &
                .AND. (z_full_c(ii,kshmin(ii)) - z_full_c(ii,ntpar_l))    &
                   <=  1.25*(zh_c(ii) - z_lcl_c(ii)) .AND.                &
                 Dt_dens_parc_Tmin(ii)  <=  0.55*Dt_dens_parc_T(ii) ) THEN

              shallow_c(ii) = .TRUE.
              ! may be problem with ntpar diagnosis for deep if wadv test sets
              ! L_shallow  false
            END IF

          ELSE
            IF ( w_avg(ii)  <  cvdiag_sh_wtest .AND.                  &
             (z_full_c(ii,ntpar_l)  <=  2500.0 .OR.                   &
                 t(ii,ntpar_l)  >=  tm)                               &
            .AND. (z_full_c(ii,kshmin(ii)) - z_full_c(ii,ntpar_l))    &
               <=  1.25*(zh_c(ii) - z_lcl_c(ii)) .AND.                &
              Dt_dens_parc_Tmin(ii)  <=  0.55*Dt_dens_parc_T(ii) ) THEN

              shallow_c(ii) = .TRUE.
              ! may be problem with ntpar diagnosis for deep if wadv test sets
              ! L_shallow  false

            END IF
          END IF  ! iconv_deep

        END IF       ! test on shmin

      END IF        ! test on cumulus

    END DO          ! ii loop

  END SELECT  ! test on type of testing inversion of not.

  !==============================================================================

ELSE    ! Some form of congestus diagnosis iconv_congestus > 1
        ! Only valid option is 1 - mass flux

  !==============================================================================

  SELECT CASE(cvdiag_inv)

  CASE (0)    ! no inversion testing
    !CDIR NODEP
    DO ii=1,nunstable
      i = index_i(ii)
      j = index_j(ii)

      IF ( cumulus_c(ii) ) THEN

        !-----------------------------------------------------------------------
        ! Conditions for congestus ;
        !    T(at top) >= freezing level ?
        !    height of top > 2.5km
        !    No need for air to be descending above ?
        !
        ! Conditions for shallow if congestus diagnosis are;
        !
        !    height of top < 2.5km
        !    Air to be descending above
        !
        ! Net result for tests mean shallow plus congestus points > shallow
        ! from old method.
        !-----------------------------------------------------------------------

        ntpar_l = ntml_c(ii)

        IF ( (t(ii,ntpar_l)  >=  tm) .AND.                                  &
             (z_full_c(ii,ntpar_l) >  2500.0)) THEN

          l_congestus(i,j) = .TRUE.

        ELSE IF (z_full_c(ii,ntpar_l) <= 2500.0) THEN

          ! Top of convection below 2500km
          IF (l_dilute) THEN   ! no w test over land
            IF (frac_land(ii) > 0.8) THEN      ! land points
              shallow_c(ii) = .TRUE.

              IF (t(ii,k_neutral(ii)) < tm) THEN
                shallow_c(ii) = .FALSE.
              END IF

            ELSE         ! ocean/coastal points do w test

              IF (w_avg(ii)  <  cvdiag_sh_wtest) THEN
                shallow_c(ii) = .TRUE.
              ELSE
                l_congestus(i,j) = .TRUE.
              END IF
            END IF

          ELSE    ! undilute ascents (icvdiag=1)

            ! Air above shallow convection less than a set value
            IF (w_avg(ii)  <  cvdiag_sh_wtest) THEN
              shallow_c(ii) = .TRUE.
            ELSE   ! Air rising therefore class as congestus rather than deep
              l_congestus(i,j) = .TRUE.
            END IF    ! w_avg test

          END IF      ! (l_dilute)
        END IF        ! test on TM etc

        !  l_congestus - all congestus points
        !  l_congestus2 - points with descending air

        IF (l_congestus(i,j)) THEN
          IF (w_avg(ii)  <  cvdiag_sh_wtest ) THEN
            l_congestus2(i,j) = .TRUE.
          END IF     ! descending air
        END IF       ! on congestus

      END IF     ! test on cumulus

    END DO      ! ii loop

  CASE (1)    ! original inversion test
    !CDIR NODEP
    DO ii=1,nunstable
      i = index_i(ii)
      j = index_j(ii)

      IF ( cumulus_c(ii) ) THEN

        !-----------------------------------------------------------------------
        ! Conditions for congestus ;
        !   As for shallow ie shmin = .TRUE.
        !    T(at top) >= freezing level ?
        !    height of top > 2.5km  < 3.5km?
        !    No need for air to be descending above ?
        !
        ! Conditions for shallow if congestus diagnosis are;
        !
        !    height of top < 2.5km
        !    Air to be descending above
        !
        ! Net result for tests mean shallow plus congestus points > shallow
        ! from old method.
        !-----------------------------------------------------------------------
        IF ( shmin(ii) ) THEN
          ntpar_l = ntml_c(ii)
          IF ( (z_full_c(ii,kshmin(ii)) - z_full_c(ii,ntpar_l))          &
                    <=  1.25*(zh_c(ii) - z_lcl_c(ii)) .AND.              &
                 Dt_dens_parc_Tmin(ii)  <=  0.55*Dt_dens_parc_T(ii) )    &
            THEN

            ! top of convection above 2500km and still below freezing level
            ! No check on descending air

            IF ( (t(ii,ntpar_l)  >=  tm) .AND.                           &
               (z_full_c(ii,ntpar_l) >  2500.0)) THEN
              l_congestus(i,j) = .TRUE.
            ELSE

              ! Top of convection below 2500km

              IF (z_full_c(ii,ntpar_l) <= 2500.0) THEN
                ! Air above shallow convection less than a set value
                IF (w_avg(ii)  <  cvdiag_sh_wtest) THEN
                  shallow_c(ii) = .TRUE.
                ELSE   ! Air rising therefore class as congestus
                       ! rather than deep
                  l_congestus(i,j) = .TRUE.
                END IF    ! w_avg test
              END IF      ! test on 2500m

            END IF        ! test on TM etc

          END IF         ! shmin

          !  l_congestus - all congestus points
          !  l_congestus2 - points with descending air

          IF (l_congestus(i,j)) THEN
            IF (w_avg(ii)  <  cvdiag_sh_wtest ) THEN
              l_congestus2(i,j) = .TRUE.
            END IF     ! descending air
          END IF       ! on congestus

        END IF       ! test on shmin

      END IF     ! test on cumulus

    END DO      ! ii loop

  END SELECT  ! test on type of testing inversion of not.

  !=======================================================================

END IF            ! test on congestus diagnosis

!=======================================================================

DO ii=1, nunstable
  i = index_i(ii)
  j = index_j(ii)
  IF (cumulus_c(ii)) THEN
    !-------------------------------------------------------------------
    ! Set mixed layer depth to z_lcl
    !-------------------------------------------------------------------
    IF (p_lcl(ii)  <   (P_theta_lev_c(ii,nlcl_c(ii)+1))) THEN
      !-------------------------------------------------------------------
      ! If LCL is diagnosed in the upper half of the layer set z_lcl to
      ! the height of the upper layer interface
      ! (in code above LCL is always set to the lower interface).
      !-------------------------------------------------------------------
      nlcl_c(ii) = nlcl_c(ii)+1
    END IF

  ELSE      ! not cumulus

    !---------------------------------------------------------------------
    !      If not cumulus, reset parameters to within bl_levels
    !---------------------------------------------------------------------
    IF (ntml_c(ii)  >   mbl) THEN
      ntml_c(ii)  = mbl
      zh_c(ii)    = z_half_c(ii,mbl+1)
    END IF

  END IF        ! test on cumulus

  !---------------------------------------------------------
  ! nlcl is not permitted to be less than nlcl_min
  !---------------------------------------------------------

  nlcl_c(ii) = MAX(nlcl_min(ii), nlcl_c(ii))

  IF ( ntml_c(ii)-nlcl_c(ii) <=2 ) THEN
    ! Cloud layer now too shallow so rediagnose as well-mixed
    cumulus_c(ii) = .FALSE.
    shallow_c(ii) = .FALSE.
    l_congestus(i,j) = .FALSE.
  END IF

END DO       ! ii loop

!----------------------------------------------------------------------
! Expand back up to full arrays
!----------------------------------------------------------------------
DO ii=1, nunstable
  i = index_i(ii)
  j = index_j(ii)

  delthvu(i,j) = delthvu_c(ii)
  cape(i,j)    = CAPE_c(ii)
  cin(i,j)     = CIN_c(ii)
  ql_ad(i,j)   = ql_ad_c(ii)

  ntpar(i,j)   = ntml_c(ii)   ! holds parcel top level
  zh(i,j)      = zh_c(ii)     ! holds parcel top at this point
  zhpar(i,j)   = zh_c(ii)

  IF (cumulus_c(ii)) THEN
    ntml(i,j)  = nlcl_c(ii)   ! now sets to LCL
    zh(i,j)    = z_half_c(ii,nlcl_c(ii)+1)    ! reset to zlcl
  ELSE
    ntml(i,j)  = ntml_c(ii)
  END IF
  nlcl(i,j)  = nlcl_c(ii)

  z_lcl_nlcl(i,j) = z_half_c(ii,nlcl_c(ii)+1)  ! LCL height
  z_lcl_uv(i,j)   = z_full_c(ii,nlcl_c(ii))    ! LCL height on uv

  ! Cumulus points always have nlcl <= MBL from previous checks
  IF (nlcl_c(ii)  >   mbl) THEN    ! only applied if not cumulus
    nlcl(i,j)    = mbl
    z_lcl_nlcl(i,j) = zh_c(ii)
    z_lcl_uv(i,j)   = z_full_c(ii,mbl-1)
  END IF

END DO       ! ii loop

!      If cumulus has been diagnosed but delthvu is negative, reset
!      cumulus and L_shallow to FALSE but leave zh and ntml at LCL

DO ii=1, nunstable
  IF (cumulus_c(ii) .AND. delthvu_c(ii)  <=  0.0) THEN
    cumulus_c(ii)    = .FALSE.
    shallow_c(ii)    = .FALSE.
  END IF
END DO

! Additional requirement for grey zone (so high res) that there is
! low level convergence (ie significant positive w) before convection
! scheme is triggered (note deep and mid-level are disabled)
IF (cldbase_opt_sh == sh_grey_closure) THEN
  DO ii=1, nunstable
    ! Start by making all convection shallow
    shallow_c(ii) = cumulus_c(ii)
    i = index_i(ii)
    j = index_j(ii)
    k = nlcl_c(ii)+1
    IF ( cumulus_c(ii) .AND. w_copy(i,j,k) < cvdiag_sh_wtest) THEN
      shallow_c(ii) = .FALSE.   ! disables convection scheme
    END IF
  END DO
END IF

! Expand back shallow and cumulus arrays
DO ii=1, nunstable
  i = index_i(ii)
  j = index_j(ii)
  cumulus(i,j)   = cumulus_c(ii)
  l_shallow(i,j) = shallow_c(ii)
  ! No points with altered entrainment go through shallow scheme.

  IF (l_shallow(i,j) .AND. entrain_coef(i,j) > 0.0) THEN
    entrain_coef(i,j) = -99.0
  END IF

END DO

!  If cumulus has been diagnosed but delthvu is negative, reset
!      congestus to FALSE but leave zh and ntml at LCL

IF (iconv_congestus > 0) THEN
  DO ii=1, nunstable
    i = index_i(ii)
    j = index_j(ii)
    IF (cumulus_c(ii) .AND. delthvu_c(ii)  <=  0.0) THEN
      l_congestus(i,j)  = .FALSE.
      l_congestus2(i,j) = .FALSE.
    END IF
  END DO       ! ii loop
END IF            ! test on congestus diagnosis

!-----------------------------------------------------------------------
! Finalise inversion thickness calculation
!-----------------------------------------------------------------------
IF (forced_cu >= on .OR. bl_res_inv == on) THEN
  DO ii=1, nunstable
    i = index_i(ii)
    j = index_j(ii)
    IF (cumulus(i,j)) THEN
      dzh(i,j)  = rmdi
    ELSE
      ! need to reset parcel top to inversion base
      ! (ie without inversion included)
      ntml(i,j) = ntml_save(ii)
      zh(i,j)   = zh_save(ii)
      IF (ntml(i,j) > mbl) THEN
        ntml(i,j) = mbl
        zh(i,j)   = z_half_c(ii,mbl+1)
      END IF
      ! restrict inversion thickness to be less than the BL depth
      dzh(i,j) = MIN( zh(i,j), zh_itop_c(ii) - zh(i,j) )
    END IF
    ! forced_cu needs parcel water content at the top of the inversion
    IF ( dzh(i,j) > 0.0 ) THEN
      zh_itop_c(ii) = zh(i,j)+dzh(i,j)
    ELSE
      ! this will be used as a lower limit on cumulus cloud base in-cloud qcl
      zh_itop_c(ii) = zh(i,j)+300.0
    END IF
    k = ntml(i,j)
    ntinv_c(ii) = k
    DO WHILE (z_half_c(ii,k+1) < zh_itop_c(ii) .AND. k < model_levels-1)
      ntinv_c(ii) = k  ! marks top level within inversion
      k=k+1
    END DO
    qcl_inv_top(i,j) = ql_parc_c(ii,ntinv_c(ii))
  END DO
END IF


!-----------------------------------------------------------------------
!     SCM  Diagnostics all versions
!-----------------------------------------------------------------------
IF (model_type == mt_single_column) THEN

  IF (l_scmdiags(scmdiag_conv)) THEN

    ! single level fields compressed
    ! Zero output
    TmpScm2d(:,:)=0.0


    ! Fill as required by expanding field
    DO ii=1,nunstable
      i=index_i(ii)
      j=index_j(ii)
      TmpScm2d(i,j) = REAL(k_plume(ii))
    END DO
    CALL scmoutput(TmpScm2d,'k_plume',                                  &
         'model level for parcel start',' ',                            &
         t_inst,d_point,default_streams,'',routinename)

    DO ii=1,nunstable
      i=index_i(ii)
      j=index_j(ii)
      TmpScm2d(i,j) = qw_plume(ii)
    END DO
    CALL scmoutput(TmpScm2d,'qw_plume',                                 &
         'initial parcel water','kg/kg',                                &
         t_inst,d_point,default_streams,'',routinename)

    DO ii=1,nunstable
      i=index_i(ii)
      j=index_j(ii)
      TmpScm2d(i,j) = sl_plume(ii)
    END DO
    CALL scmoutput(TmpScm2d,'sl_plume',                                 &
         'initial parcel energy','J',                                   &
         t_inst,d_point,default_streams,'',routinename)

    ! model level output - compressed fields

    DO k= 1,model_levels
      DO ii=1,nunstable
        i=index_i(ii)
        j=index_j(ii)
        TmpScm3d(i,j,k)= buoyancy(ii,k)
      END DO
    END DO
    CALL scmoutput(TmpScm3d,'buoy_undil',                               &
         'parcel buoyancy -undilute','K',                               &
         t_inst,d_wet,default_streams,'',routinename)

    DO k= 1,model_levels
      DO ii=1,nunstable
        i=index_i(ii)
        j=index_j(ii)
        TmpScm3d(i,j,k)= t_parc(ii,k)
      END DO
    END DO
    CALL scmoutput(TmpScm3d,'t_parc_undil',                             &
         'Undilute Parcel temperature','K',                             &
         t_inst,d_wet,default_streams,'',routinename)

    IF (icvdiag >= 2) THEN          ! dilute parcel options only
      DO k= 1,model_levels
        DO ii=1,nunstable
          i=index_i(ii)
          j=index_j(ii)
          TmpScm3d(i,j,k)=buoyancy_dil(ii,k)
        END DO
      END DO
      CALL scmoutput(TmpScm3d,'buoy_dil',                               &
           'parcel buoyancy - dilute','K',                              &
           t_inst,d_wet,default_streams,'',routinename)

      DO k= 1,model_levels
        DO ii=1,nunstable
          i=index_i(ii)
          j=index_j(ii)
          TmpScm3d(i,j,k)=t_parc_dil(ii,k)
        END DO
      END DO
      CALL scmoutput(TmpScm3d,'t_parc_dil',                             &
           'dilute parcel temperature','K',                             &
           t_inst,d_wet,default_streams,'',routinename)

      DO k=1, model_levels
        DO ii=1, nunstable
          i=index_i(ii)
          j=index_j(ii)
          TmpScm3d(i,j,k)=entrain_fraction(ii,k)
        END DO
      END DO
      CALL scmoutput(TmpScm3d,'entrain_frac',                           &
           'entrainment fraction',' ',                                  &
           t_inst,d_wet,default_streams,'',routinename)

    END IF ! test on icvdiag > 2 i.e. dilute parcel

  END IF ! scmdiag_conv

  IF (l_scmdiags(scmdiag_bl)) THEN

    DO k=1, model_levels
      DO ii=1, nunstable
        i=index_i(ii)
        j=index_j(ii)
        TmpScm3d(i,j,k)=par_svl(ii,k)
      END DO
    END DO
    CALL scmoutput(TmpScm3d,'thv_par',                                  &
         'thetav parcel','K',                                           &
         t_inst,d_wet,default_streams,'',routinename)

    DO k=1, model_levels
      DO ii=1, nunstable
        i=index_i(ii)
        j=index_j(ii)
        TmpScm3d(i,j,k)=env_svl(ii,k)
      END DO
    END DO
    CALL scmoutput(TmpScm3d,'thv_env',                                  &
         'Environment thetav','K',                                      &
         t_inst,d_wet,default_streams,'',routinename)

  END IF ! scmdiag_bl

END IF ! model_type
!-----------------------------------------------------------------------

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE conv_diag_comp_5a

END MODULE conv_diag_comp_5a_mod
