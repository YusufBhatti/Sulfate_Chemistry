! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!  Calculate and output some SCM diagnostics

SUBROUTINE dgnstcs_scm_main                                                   &
  ( row_length, rows, land_points, tr_levels                                  &
  , tr_vars, sm_levels, nsmax                                                 &
  , st_levels, ntype, rho, timestep, u, v, t                                  &
  , theta, q, qcl, qcf, layer_cloud, cca, ccw, t_deep_soil, p_star, tstar     &
  , smc, canopy_gb, snodep, zh, z0msea, smcl, sthu, sthf, gs                  &
  , nsnow, snowl_ds, snowl_ice, snowl_liq, snowl_tmpr, snowl_rho              &
  , snowl_rgrain, snows_rgrain                                                &
  , lw_incs, photosynth_act_rad, tstar_tile, aerosol, free_tracers            &
  , p_theta_levels, p                                                         &
  , iccb, icct, w, w_adv, area_cloud_fraction, bulk_cloud_fraction            &
  , cloud_fraction_liquid, cloud_fraction_frozen, cclwp, nSCMDpkgs            &
  , L_SCMDiags )

USE cv_run_mod, ONLY: l_ccrad, l_3d_cca

USE level_heights_mod, ONLY: r_theta_levels, r_rho_levels

USE nlsizes_namelist_mod, ONLY: model_levels

USE s_scmop_mod, ONLY: stream                                               &
  , t_inst, t_avg, t_max, t_min, t_acc, t_div, t_mult, t_acc_div            &
  , t_acc_mult, t_const, only_radsteps                                      &
  , d_sl, d_soilt, d_bl, d_wet, d_all, d_soilm, d_tile, d_vis, d_point      &
  , d_allxtra, d_land, d_cloud, d_mlsnow, default_streams                   &
  , SCMDiag_gen, SCMDiag_rad, SCMDiag_bl, SCMDiag_surf, SCMDiag_land        &
  , SCMDiag_sea, SCMDiag_lsp, SCMDiag_conv, SCMDiag_lscld, SCMDiag_pc2      &
  , SCMDiag_forc, SCMDiag_incs, SCMDiag_gwd, SCMDiag_MLSnow
USE scmoutput_mod, ONLY: scmoutput

USE qsat_mod, ONLY: qsat, qsat_wat

USE parkind1, ONLY: jpim, jprb
USE yomhook,  ONLY: lhook, dr_hook

IMPLICIT NONE

! Description:
!   Essentially just a lot of calls to SCMoutput

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Single Column Model

! Code Description:
!   Language: Fortran 90

! All parameters are Intent In

! Variables passed in through the parameter list...
INTEGER ::      &
  row_length    &
, rows          &
, land_points   &
, tr_levels     &
, tr_vars       &
, sm_levels     &
, st_levels     &
, ntype         &
, nsmax


REAL ::                                   &
  rho(row_length,rows,model_levels)       &! Density*r^2 (kg/m)
, u(row_length,rows,model_levels)         &! Zonal wind (m/s)
, v(row_length,rows,model_levels)         &! Meridional wind (m/s)
, w(row_length,rows,0:model_levels)       &! Vertical velocity (m/s)
, w_adv(row_length,rows,0:model_levels)   &! Advective vertical velocity
                                           ! (m/s)
, t(row_length,rows,model_levels)         &! Temperature(K)
, theta(row_length,rows,model_levels)     &! Potential temperature (K)
, q(row_length,rows,model_levels)         &! Specific humidity (kg/kg)
, qcf(row_length,rows,model_levels)       &! Cloud ice content (kg/kg)
, qcl(row_length,rows,model_levels)       &! Cloud water content(kg/kg)
, layer_cloud(row_length,rows,model_levels) &! Layer cloud amount
                                             ! (decimal fraction)
, cca(row_length,rows,model_levels)       &! Convective cloud amount
, ccw(row_length,rows,model_levels)        ! Convective Cloud Water passed
                                           ! to radiation scheme (kg/kg)

REAL ::                                             &
  area_cloud_fraction(row_length,rows,model_levels) &! Area cloud amount
                                                     ! (decimal fraction)
, bulk_cloud_fraction(row_length,rows,model_levels) &! Cloud amount
                                                     ! (decimal fraction)
, cloud_fraction_liquid(row_length,rows,model_levels) &! Liquid cloud amount
                                                     ! (decimal fraction)
, cloud_fraction_frozen(row_length,rows,model_levels) &! Frozen cloud amount
                                                     ! (decimal fraction)

, cclwp(row_length,rows)             &! condensed water path (kg/m^2)
, t_deep_soil(land_points,st_levels) &! Deep soil temperatures (K)
, p_star(row_length,rows)            &! Pressure at earth's surface
, tstar(row_length,rows)             &! Surface temperature (K)
, smc(land_points)                   &! Soil moisture content(kg/m^2)
, smcl(land_points,sm_levels)        &! Soil moisture content in layers
                                      ! (kg/m^2)
, sthf(land_points,sm_levels)        &! Frozen soil moisture content of
                                      ! each layer as a fraction of
                                      ! saturation.
, canopy_gb(land_points)             &! Canopy water content (kg/m^2)
, snodep(row_length,rows)            &! Snow depth (kg/m^2)
, p(row_length,rows,model_levels+1)   ! Pressure on rho levels

REAL ::                                         &
  p_theta_levels(row_length,rows,model_levels)  &! Pressure on theta levels
, zh(row_length,rows)                           &! Height above surface of top
                                                 ! of boundary layer (m)
, z0msea(row_length,rows)            &! Sea surface roughness length
, sthu(land_points,sm_levels)        &! Unfrozen soil moisture content
                                      ! of each layer as a fraction of
                                      ! saturation. (kg/m^2)
, gs(row_length*rows)                &! Stomatal conductance
, LW_incs(row_length,rows,0:model_levels) &
, photosynth_act_rad(row_length,rows)     &! Net downward shortwave
                                           ! radiation in band 1 (w/m^2).
, tstar_tile(land_points,ntype)       &! Surface tile temperature
, aerosol(row_length,rows,model_levels)   &! Aerosol values ; only used if
                                           ! l_murk=.true. ; default .false.
, free_tracers(row_length,rows,tr_levels,tr_vars)    &
, timestep

! Multilayer snow diagnostics
REAL ::                     &
  nsnow(land_points,ntype)
!                                            ! Number of layers in the
!                                            ! snowpack
REAL ::                     &
  snowl_ds(land_points,ntype,nsmax)                    &
!                                            ! Thickness of each snow layer
  , snowl_ice(land_points,ntype,nsmax)                 &
!                                            ! Frozen mass of each snow layer
  , snowl_liq(land_points,ntype,nsmax)                 &
!                                            ! Unforzen mass of each snow layer
  , snowl_tmpr(land_points,ntype,nsmax)                &
!                                            ! Temperature of each snow layer
  , snowl_rho(land_points,ntype,nsmax)                 &
!                                            ! Density of each snow layer
  , snowl_rgrain(land_points,ntype,nsmax)
!                                            ! Grain sice in each snow layer
REAL ::                     &
  snows_rgrain(land_points,ntype)
!                                            ! Grain sice on surface
    
INTEGER ::                  &
  iccb(row_length, rows)    &! Convective cloud base and top
, icct(row_length, rows)     ! at levels 1 to model_levels

INTEGER :: nSCMDpkgs      ! No of SCM diagnostics packages
LOGICAL :: L_SCMDiags(nSCMDpkgs)  ! Diagnostics packages logicals

! Local variables...

CHARACTER(LEN=*), PARAMETER ::  RoutineName = 'DGNSTCS_SCM_MAIN'
CHARACTER (LEN=2)        :: cdiv    ! Character string for tracer counter

REAL ::                            &
  sum_p_col(row_length,rows)       &! Sums of pressures on theta levels
, coltemp(row_length,rows)         &! Average pressure-weighted
                                    ! column temperature
, colq(row_length,rows)            &! Average pressure-weighted
                                    ! column humidity
, accum_pptn(row_length,rows)      &! Accumulated precipitation
, pptn_rate(row_length,rows)       &! Precipitation rate
, ccb_m(row_length,rows)           &! Convective cloud base (m)
, ccb_Pa(row_length,rows)          &! Convective cloud base (Pa)
, cct_m(row_length,rows)           &! Convective cloud top  (m)
, cct_Pa(row_length,rows)          &! Convective cloud top  (Pa)
, rh(row_length,rows,model_levels,2) &! Relative humidity (over land and water)
, qst                               ! A saturation mixing ratio


REAL ::                                  &
  z_theta(row_length,rows,model_levels)  &! height above surface (m) (Theta)
, z_rho(row_length,rows,model_levels)    &! height above surface (m) (Rho)
, rho_only(row_length,rows,model_levels)  ! Actual density (kg/m3)

REAL :: a2out(row_length,rows,model_levels)
INTEGER :: i,j,k,n,l

! Dr Hook
!==============================
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb) :: zhook_handle

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)

! Store diagnostics...

!-----------------------------------------------------------------------
!     SCM General Diagnostics Package
!-----------------------------------------------------------------------
IF (L_SCMDiags(SCMDiag_gen)) THEN

  ! Stash 0,002
  CALL scmoutput(u, 'u'                                                     &
    , 'Zonal wind', 'm/s'                                                   &
    , t_avg, d_all, default_streams, '', RoutineName)

  ! Stash 0,003
  CALL scmoutput(v, 'v'                                                     &
    , 'Meridional wind', 'm/s'                                              &
    , t_avg, d_all, default_streams, '', RoutineName)

  ! Stash 30,111
  CALL scmoutput(t, 'T'                                                     &
    , 'Temperature', 'K'                                                    &
    , t_avg, d_all, default_streams, '', RoutineName)

  ! Stash 0,004
  CALL scmoutput(theta, 'theta'                                             &
    , 'Potential temperature', 'K'                                          &
    , t_avg, d_all, default_streams, '', RoutineName)

  ! Stash 0,010
  CALL scmoutput(q, 'q'                                                     &
    , 'Specific humidity', 'kg/kg'                                          &
    , t_avg, d_wet, default_streams, '', RoutineName)

  DO k=1,model_levels
    DO j=1,rows
      DO i=1,row_length
        a2out(i,j,k) = w(i,j,k)
      END DO
    END DO
  END DO

  ! Stash 0,150
  CALL scmoutput(a2out, 'w'                                                 &
    , 'Vertical velocity', 'm/s'                                            &
    , t_avg, d_all, default_streams, '', RoutineName)

  ! Stash 0,407
  CALL scmoutput(p, 'p_rho'                                                 &
    , 'Pressure on rho levels', 'Pa', t_avg, d_all, default_streams, ''     &
    , RoutineName)

  ! Stash 0,408
  CALL scmoutput(p_theta_levels, 'p_theta'                                  &
    , 'Pressure on theta levels', 'Pa'                                      &
    , t_avg, d_all, default_streams, '', RoutineName)

  DO k=1, model_levels
    DO j=1, rows
      DO i=1, row_length
        ! Height above model surface as calculation in ni_conv_ctl
        z_theta(i,j,k)  = r_theta_levels(i,j,k) - r_theta_levels(i,j,0)
        z_rho(i,j,k)    = r_rho_levels(i,j,k)   - r_theta_levels(i,j,0)
        rho_only(i,j,k) = rho(i,j,k)                                        &
                        / (r_rho_levels(i,j,k) * r_rho_levels(i,j,k))
      END DO ! i
    END DO ! j
  END DO ! k

  ! Stash 15,101
  CALL scmoutput(z_theta, 'h_theta'                                         &
    , 'Height of model theta levels', 'm'                                   &
    , t_avg, d_all, default_streams, '', RoutineName)

  ! Stash 15,102
  CALL scmoutput(z_rho, 'h_rho'                                             &
    , 'Height of model rho levels', 'm'                                     &
    , t_avg, d_all, default_streams, '', RoutineName)

  ! Stash 0,253
  CALL scmoutput(rho, 'rho_r2'                                              &
    , 'Density *r*r after timestep', 'kg/m'                                 &
    , t_avg, d_all, default_streams, '', RoutineName)

  CALL scmoutput(rho_only, 'rho_only'                                       &
    , 'Density after timestep', 'kg/m3'                                     &
    , t_avg, d_all, default_streams, '', RoutineName)

  ! Calculate relative humidities on theta levels
  DO j=1, rows
    DO i=1, row_length
      DO k=1, model_levels

        CALL qsat(qst, t(i,j,k), p_theta_levels(i,j,k))
        rh(i,j,k,1) = q(i,j,k) / qst

        CALL qsat_wat(qst, t(i,j,k), p_theta_levels(i,j,k))
        rh(i,j,k,2) = q(i,j,k) / qst

      END DO ! k
    END DO ! i
  END DO ! j

  ! Stash 30,113
  CALL scmoutput(rh(1,1,1,1), 'rh'                                          &
    , 'Relative humidity', '%'                                              &
    , t_avg, d_wet, default_streams, '', RoutineName)

  CALL scmoutput(rh(1,1,1,2), 'rh2'                                         &
    , 'Relative humidity over liquid water', '%'                            &
    , t_avg, d_wet, default_streams, '', RoutineName)

  ! Calculate pressure-weighted average temperature and humidity
  i = 1
  j = 1
  colq(i,j) = 0.0
  coltemp(i,j) = 0.0
  sum_p_col(i,j) = 0.0
  DO k=1, model_levels
    sum_p_col(i,j) = sum_p_col(i,j) + p_theta_levels(i,j,k)
    coltemp(i,j)   = coltemp(i,j)   + t(i,j,k)*p_theta_levels(i,j,k)
    colq(i,j)      = colq(i,j)      + q(i,j,k)*p_theta_levels(i,j,k)
  END DO ! k

  CALL scmoutput(sum_p_col, 'sum_p_col'                                     &
    , 'Sum of theta-level pressures', 'Pa'                                  &
    , t_acc, d_sl, default_streams, '', RoutineName)

  CALL scmoutput(coltemp, 'tatmos'                                          &
    , 'Mean column temperature', 'K'                                        &
    , t_acc_div, d_sl, default_streams, 'sum_p_col', RoutineName)

  CALL scmoutput(colq, 'qatmos'                                             &
    , 'Mean column specific humidity', 'kg/kg'                              &
    , t_acc_div, d_sl, default_streams, 'sum_p_col', RoutineName)

  CALL scmoutput(aerosol, 'aerosol'                                         &
    , 'Murk aerosol concentration', 'micro-g/kg'                            &
    , t_avg, d_wet, default_streams, '', RoutineName)

  IF (tr_vars > 0) THEN

    DO n=1,tr_vars

      a2out(:,:,:) = 0.0
      DO k=1, tr_levels
        a2out(:,:,k) = free_tracers(:,:,k,n)
      END DO
      WRITE(cdiv,'(I2.2)') n

      CALL SCMoutput(a2out,                                                 &
         'tracer'//cdiv,'Concentration of tracer '//cdiv,'kg/kg',           &
         t_avg,d_all,default_streams,'',RoutineName)

    END DO

  END IF  ! test on tr_vars

END IF ! L_SCMDiags(SCMDiag_gen)

!-----------------------------------------------------------------------
!     SCM General OR Large Scale Cloud Diagnostics Package
!-----------------------------------------------------------------------
IF (L_SCMDiags(SCMDiag_gen) .OR. L_SCMDiags(SCMDiag_lscld)) THEN

  ! Stash 0,254
  CALL scmoutput(qcl, 'qcl'                                                 &
    , 'QCL cloud water content', 'kg/kg'                                    &
    , t_avg, d_wet, default_streams, '', RoutineName)

  ! Stash 0,012
  CALL scmoutput(qcf, 'qcf'                                                 &
    , 'QCF cloud ice content', 'kg/kg'                                      &
    , t_avg, d_wet, default_streams, '', RoutineName)

  CALL scmoutput(layer_cloud, 'layer_cloud'                                 &
    , 'Layer cloud amount', 'Fraction'                                      &
    , t_avg, d_wet, default_streams, '', RoutineName)

  ! Stash 0, 265
  CALL scmoutput(area_cloud_fraction, 'acf'                                 &
    , 'Area cloud fraction', 'Fraction'                                     &
    , t_avg, d_wet, default_streams, '', RoutineName)

  ! Stash 0, 266
  CALL scmoutput(bulk_cloud_fraction, 'bcf'                                 &
    , 'Bulk cloud fraction', 'Fraction'                                     &
    , t_avg, d_wet, default_streams, '', RoutineName)

  ! Stash 0, 267
  CALL scmoutput(cloud_fraction_liquid, 'cfl'                               &
    , 'Liquid cloud fraction', 'Fraction'                                   &
    , t_avg, d_wet, default_streams, '', RoutineName)

  ! Stash 0, 268
  CALL scmoutput(cloud_fraction_frozen, 'cff'                               &
    , 'Frozen cloud fraction', 'Fraction'                                   &
    , t_avg, d_wet, default_streams, '', RoutineName)

END IF ! L_SCMDiags(SCMDiag_gen) .OR. L_SCMDiags(SCMDiag_lscld)

!-----------------------------------------------------------------------
!     SCM General OR Convection Diagnostics Package
!-----------------------------------------------------------------------
IF (L_SCMDiags(SCMDiag_gen) .OR. L_SCMDiags(SCMDiag_conv)) THEN

  ! Stash 0,016
  ! Condensed water path convective cloud passed to radiation
  CALL scmoutput(cclwp, 'cclwp2rad'                                         &
    , 'cclwp passed to Radiation', 'kg/m^2'                                 &
    , t_avg, d_sl, default_streams, '', RoutineName)


  IF (l_3d_cca) THEN
    ! Stash 0,211
    ! Convective Cloud Amount passed to Radiation
    CALL SCMoutput(cca                                                      &
      , 'cca2rad','CCA passed to radiation', 'Fraction'                     &
      , t_avg, d_all, default_streams, '', RoutineName)
  ELSE
    ! Stash 0,211
    ! Convective Cloud Amount passed to Radiation
    CALL SCMoutput(cca(1,1,1)                                               &
      , 'cca2rad','CCA passed to radiation', 'Fraction'                     &
      , t_avg, d_sl, default_streams, '', RoutineName)
  END IF

  ! Stash 0,212
  ! Convective Cloud Water passed to Radiation
  IF (l_ccrad) THEN

    CALL scmoutput(ccw, 'ccw2rad'                                           &
      , 'CCW passed to radiation', 'kg/kg'                                  &
      , t_avg, d_wet, default_streams, '', RoutineName)
  END IF

  DO k=1,model_levels
    DO j=1,rows
      DO i=1,row_length
        a2out(i,j,k) = w_adv(i,j,k)
      END DO
    END DO
  END DO

  ! Stash 0,258
  CALL scmoutput(a2out, 'w_adv'                                             &
    , 'Advective vertical velocity', 'm/s'                                  &
    , t_avg, d_all, default_streams, '', RoutineName)

END IF ! L_SCMDiags(SCMDiag_gen) .OR. L_SCMDiags(SCMDiag_conv)

!-----------------------------------------------------------------------
!     SCM General OR Surface Diagnostics Package
!-----------------------------------------------------------------------
IF (L_SCMDiags(SCMDiag_gen) .OR. L_SCMDiags(SCMDiag_surf)) THEN

  ! Stash 0,409
  CALL scmoutput(p_star, 'pstar'                                            &
    , 'Surface pressure', 'Pa'                                              &
    , t_avg, d_sl, default_streams, '', RoutineName)

  ! Stash 0,024
  CALL scmoutput(tstar, 'tstar'                                             &
    , 'Surface temperature', 'K'                                            &
    , t_avg, d_sl, default_streams, '', RoutineName)

END IF ! L_SCMDiags(SCMDiag_gen) .OR. L_SCMDiags(SCMDiag_surf)

!-----------------------------------------------------------------------
!     SCM General OR Radiation Diagnostics Package
!-----------------------------------------------------------------------
IF (L_SCMDiags(SCMDiag_gen) .OR. L_SCMDiags(SCMDiag_rad)) THEN

  ! Stash 2,232
  CALL scmoutput(lw_incs, 'lwrate_day'                                      &
    , 'LW heating rate', 'K/day'                                            &
    , t_mult, d_allxtra, default_streams, 'ntspday', RoutineName)

  CALL scmoutput(lw_incs, 'lwrate_step'                                     &
    , 'LW heating rate', 'K/timestep'                                       &
    , t_avg, d_allxtra, default_streams, 'ntspday', RoutineName)

  CALL scmoutput(lw_incs, 'lwrate_acc'                                      &
    , 'Accumulated LW heating rate over dumping period', 'K/ts*dumps'       &
    , t_acc, d_allxtra, default_streams, '', RoutineName)

END IF ! L_SCMDiags(SCMDiag_gen) .OR. L_SCMDiags(SCMDiag_rad)

!-----------------------------------------------------------------------
!     SCM Surface OR Radiation Diagnostics Package
!-----------------------------------------------------------------------
IF (L_SCMDiags(SCMDiag_surf) .OR. L_SCMDiags(SCMDiag_rad)) THEN

  CALL scmoutput(photosynth_act_rad, 'down_surf_sw_b1'                      &
    , 'Downward SW radn in band 1', 'W/m^2'                                 &
    , t_avg, d_sl, default_streams, '', RoutineName)

END IF ! L_SCMDiags(SCMDiag_surf) .OR. L_SCMDiags(SCMDiag_rad)

!-----------------------------------------------------------------------
!     SCM Convection Diagnostics Package
!-----------------------------------------------------------------------
IF (L_SCMDiags(SCMDiag_conv)) THEN

  ! Convert ccb as a model level to ccb in meters and Pascals
  DO j=1, rows
    DO i=1, row_length
      IF (iccb(i,j) /= 0) THEN
        ccb_m(i,j)  = r_theta_levels(i,j,iccb(i,j)) - r_theta_levels(i,j,0)
        cct_m(i,j)  = r_theta_levels(i,j,icct(i,j)) - r_theta_levels(i,j,0)
        ccb_pa(i,j) = p_theta_levels(i,j,iccb(i,j))
        cct_pa(i,j) = p_theta_levels(i,j,icct(i,j))
      ELSE
        ccb_m(i,j)  = 0.0
        cct_m(i,j)  = 0.0
        ccb_pa(i,j) = 0.0
        cct_pa(i,j) = 0.0
      END IF
    END DO ! i
  END DO ! j

  CALL scmoutput                                                            &
     ( ccb_m*cca(1:row_length,1:rows, 1), 'ccb_z'                           &
     , 'Convective cloud base height (weighted average)', 'm'               &
     , t_div, d_sl, default_streams, 'cca2rad', RoutineName )

  CALL scmoutput                                                            &
     ( ccb_pa*cca(1:row_length,1:rows, 1), 'ccb_pa'                         &
     , 'Convective cloud base pressure (weighted average)', 'Pa'            &
     , t_div, d_sl, default_streams, 'cca2rad', RoutineName )

  CALL scmoutput                                                            &
     ( cct_m*cca(1:row_length,1:rows,1), 'cct_z'                            &
     , 'Convective cloud top height (weighted average)', 'm'                &
     , t_div, d_sl, default_streams, 'cca2rad', RoutineName )

  CALL scmoutput                                                            &
     ( cct_pa*cca(1:row_length,1:rows,1), 'cct_pa'                          &
     , 'Convective cloud top pressure (weighted average)', 'Pa'             &
     , t_div, d_sl, default_streams, 'cca2rad', RoutineName )

END IF ! L_SCMDiags(SCMDiag_conv)

!-----------------------------------------------------------------------
!     SCM Boundary Layer Diagnostics Package
!-----------------------------------------------------------------------
IF (L_SCMDiags(SCMDiag_bl)) THEN

  CALL scmoutput(zh, 'bl_depth'                                             &
    , 'Boundary layer depth (zh)', 'm'                                      &
    , t_avg, d_sl, default_streams, '', RoutineName)

END IF ! L_SCMDiags(SCMDiag_bl)

!-----------------------------------------------------------------------
!     SCM Sea Points Diagnostics Package
!-----------------------------------------------------------------------
IF (L_SCMDiags(SCMDiag_sea)) THEN

  CALL scmoutput(z0msea, 'sea_roughness'                                    &
    , 'Sea surface roughness length', 'm'                                   &
    , t_avg, d_sl, default_streams, '', RoutineName)

END IF ! L_SCMDiags(SCMDiag_sea)

!-----------------------------------------------------------------------
!     SCM Land Points Diagnostics Package
!-----------------------------------------------------------------------

IF (L_SCMDiags(SCMDiag_land) .AND. land_points > 0) THEN

  ! Stash 0,020
  CALL scmoutput(t_deep_soil, 'tsoildeep'                                   &
    , 'Deep soil temperatures', 'K'                                         &
    , t_avg, d_soilt, default_streams, '', RoutineName)

  ! Stash 8,208
  CALL scmoutput(smc, 'soilmoist'                                           &
    , 'Soil moisture content', 'kg/m^2'                                     &
    , t_avg, d_sl, default_streams, '', RoutineName)

  ! Stash 8,223
  CALL scmoutput(smcl, 'soilmoistlay'                                       &
    , 'Layer soil moisture content', 'kg/m^2'                               &
    , t_avg, d_soilm, default_streams, '', RoutineName)

  ! Stash 0,214
  CALL scmoutput(sthu, 'soilmoistunfroz'                                    &
    , 'Unfrozen soil moisture content', 'kg/m^2'                            &
    , t_avg, d_soilm, default_streams, '', RoutineName)

  ! Stash 0,215
  CALL scmoutput(sthf, 'soilmoistfroz'                                      &
    , 'Frozen soil moisture content', 'kg/m^2'                              &
    , t_avg, d_soilm, default_streams, '', RoutineName)

  CALL scmoutput(gs, 'stoma_cond'                                           &
    , 'Stomatal conductance', 'm/s'                                         &
    , t_avg, d_sl, default_streams, '', RoutineName)

  ! Stash 0,022
  CALL scmoutput(canopy_gb, 'canopy_gb'                                     &
    , 'Canopy water content', 'kg/m^2'                                      &
    , t_avg, d_sl, default_streams, '', RoutineName)

  ! Stash 0,023
  CALL scmoutput(snodep, 'snowdepth'                                        &
    , 'Snow depth', 'kg/m^2'                                                &
    , t_avg, d_sl, default_streams, '', RoutineName)

  ! Stash 3,316
  CALL scmoutput(tstar_tile, 'tstartl'                                      &
    , 'Tile surface temp', 'K'                                              &
    , t_avg, d_tile, default_streams, '', RoutineName)

END IF ! L_SCMDiags(SCMDiag_land)/land_points

! ----------------------------------------------------------------------
! Multilayer snow diagnostics
! ----------------------------------------------------------------------
IF (L_SCMDiags(SCMDiag_MLSnow) .AND. land_points > 0 .AND.                  &
     nsmax > 0 ) THEN

  ! Stash 0,231
  CALL scmoutput(snows_rgrain, 'snows_rgrain'                               &
      , 'Surface snow grain size on tiles', 'Microns'                       &
      , t_avg, d_tile, default_streams, '', RoutineName)

  ! Stash 0,380
  CALL scmoutput(nsnow, 'nsnow'                                             &
      , 'Number of snow layers on tiles', 'None'                            &
      , t_avg, d_tile, default_streams, '', RoutineName)

  ! Stash 0,381
  CALL scmoutput(snowl_ds, 'snowl_ds'                                     &
    , 'Thicknesses of snow layers on tiles', 'm'                          &
    , t_avg, d_mlsnow, default_streams, '', RoutineName)

  ! Stash 0,382
  CALL scmoutput(snowl_ice, 'snowl_ice'                                   &
    , 'Frozen mass of snow layers on tiles', 'kg.m-2'                     &
    , t_avg, d_mlsnow, default_streams, '', RoutineName)

  ! Stash 0,383
  CALL scmoutput(snowl_liq, 'snowl_liq'                                   &
    , 'Unfrozen mass of snow layers on tiles', 'kg.m-2'                   &
    , t_avg, d_mlsnow, default_streams, '', RoutineName)

  ! Stash 0,384
  CALL scmoutput(snowl_tmpr, 'snowl_tmpr'                                 &
    , 'Temperature of snow layers on tiles', 'K'                          &
    , t_avg, d_mlsnow, default_streams, '', RoutineName)

  ! Stash 0,385
  CALL scmoutput(snowl_rho, 'snowl_rho'                                   &
    , 'Density of snow layers on tiles', 'kg.m-3'                         &
    , t_avg, d_mlsnow, default_streams, '', RoutineName)

  ! Stash 0,386
  CALL scmoutput(snowl_rgrain, 'snowl_rgrain'                             &
    , 'Grain size in snow layers on tiles', 'um'                          &
    , t_avg, d_mlsnow, default_streams, '', RoutineName)

END IF ! L_SCMDiags(SCMDiag_MLSnow)/land_points/nsmax

IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)

RETURN
END SUBROUTINE dgnstcs_scm_main

